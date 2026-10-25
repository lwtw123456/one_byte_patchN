#ifndef ONE_BYTE_PATCHN_HPP
#define ONE_BYTE_PATCHN_HPP

#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <span>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace one_byte_patchN {

struct result {
    bool found = false;
    std::size_t byte_index = 0;
    unsigned char new_byte = 0;
    mpz_class new_n;
    mpz_class new_d;
    mpz_class e;
    double seconds = 0.0;
};

namespace detail {

inline long long now_ns() {
    return static_cast<long long>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
}

inline std::uint32_t sieve_limit_for(unsigned int nbits) {
    return nbits <= 2048 ? 30000u : 300000u;
}

inline std::vector<unsigned char> to_bytes_be(const mpz_class& v) {
    std::vector<unsigned char> out;
    const std::size_t nbytes = (mpz_sizeinbase(v.get_mpz_t(), 2) + 7) / 8;
    out.resize(nbytes);
    std::size_t count = 0;
    mpz_export(out.data(), &count, 1, 1, 1, 0, v.get_mpz_t());
    out.resize(count);
    return out;
}

template <std::size_t Limit>
struct prime_table {
    static constexpr std::size_t cap_big = Limit / 8 + 64;
    static constexpr std::size_t cap_small = 64;
    std::array<std::uint32_t, cap_big> big{};
    std::array<std::uint8_t, cap_big> big_inv{};
    std::size_t size_big = 0;
    std::array<std::uint32_t, cap_small> small{};
    std::array<std::uint8_t, cap_small> small_inv{};
    std::size_t size_small = 0;
};

template <std::size_t Cap>
constexpr prime_table<Cap> build_prime_table(std::size_t limit) {
    std::array<std::uint64_t, (Cap / 2) / 64 + 2> comp{};
    constexpr auto get = [](const auto& c, std::size_t i) {
        return (c[i >> 6] >> (i & 63)) & 1ull;
    };
    for (std::size_t i = 3; i * i <= limit; i += 2) {
        if (get(comp, i / 2)) continue;
        for (std::size_t j = i * i; j <= limit; j += 2 * i)
            comp[(j / 2) >> 6] |= 1ull << ((j / 2) & 63);
    }

    prime_table<Cap> t{};
    for (std::size_t i = 3; i <= limit; i += 2) {
        if (get(comp, i / 2)) continue;
        const std::uint32_t p = static_cast<std::uint32_t>(i);
        std::uint8_t inv = 1;
        for (int k = 0; k < 3; ++k)
            inv = static_cast<std::uint8_t>(inv * (2u - p * inv));
        if (p > 255) {
            t.big[t.size_big] = p;
            t.big_inv[t.size_big] = inv;
            ++t.size_big;
        } else {
            t.small[t.size_small] = p;
            t.small_inv[t.size_small] = inv;
            ++t.size_small;
        }
    }
    return t;
}

// 0 = 运行期建表（exe 小 206 KB，首次 find() 多 0.05~0.6 ms）
// 1 = 编译期建表（表进 .rdata，零建表成本）
#ifndef BYTE_PATCH_TABLE_AT_COMPILE_TIME
#define BYTE_PATCH_TABLE_AT_COMPILE_TIME 0
#endif

#if BYTE_PATCH_TABLE_AT_COMPILE_TIME
inline constexpr auto kPrimes30k = build_prime_table<30000>(30000);
inline constexpr auto kPrimes300k = build_prime_table<300000>(300000);
#else
// volatile 读：防止建表被折叠回 .rdata
template <std::size_t Cap>
std::size_t runtime_table_limit() {
    static volatile std::size_t seed = Cap;
    return seed;
}

// 堆上建表：函数内 static 会被放进 .data，仍占 exe 体积
template <std::size_t Cap>
const prime_table<Cap>& prime_table_instance() {
    static const auto* t = new prime_table<Cap>(build_prime_table<Cap>(runtime_table_limit<Cap>()));
    return *t;
}
#endif

struct prime_view {
    std::span<const std::uint32_t> big;
    std::span<const std::uint8_t> big_inv;
    std::span<const std::uint32_t> small;
    std::span<const std::uint8_t> small_inv;
};

template <class Table>
inline prime_view view_of(const Table& t) {
    return {std::span<const std::uint32_t>(t.big.data(), t.size_big),
            std::span<const std::uint8_t>(t.big_inv.data(), t.size_big),
            std::span<const std::uint32_t>(t.small.data(), t.size_small),
            std::span<const std::uint8_t>(t.small_inv.data(), t.size_small)};
}

inline prime_view primes_for(std::uint32_t limit) {
#if BYTE_PATCH_TABLE_AT_COMPILE_TIME
    return limit <= 30000 ? view_of(kPrimes30k) : view_of(kPrimes300k);
#else
    return limit <= 30000 ? view_of(prime_table_instance<30000>())
                          : view_of(prime_table_instance<300000>());
#endif
}

class sieve {
public:
    void init(const mpz_class& n, std::uint32_t limit, bool n_odd) {
        pv_ = primes_for(limit);
        const std::size_t nb = pv_.big.size();
        t_big_.resize(nb);
        for (std::size_t i = 0; i < nb; ++i)
            t_big_[i] = static_cast<std::uint32_t>(mpz_fdiv_ui(n.get_mpz_t(), pv_.big[i]));

        const std::size_t ns = pv_.small.size();
        for (std::size_t i = 0; i < ns; ++i)
            t_small_[i] = static_cast<std::uint32_t>(mpz_fdiv_ui(n.get_mpz_t(), pv_.small[i]));

        n_odd_ = n_odd;
    }

    std::size_t prime_count() const { return pv_.big.size() + pv_.small.size(); }

    void scan_position(std::size_t m, std::uint8_t old_b,
                       std::uint8_t cnt[256], std::uint32_t killer[256]) {
        std::memset(cnt, 0, 256);
        std::memset(killer, 0, 256 * sizeof(std::uint32_t));

        auto mark = [&](int nb, std::uint32_t p) {
            if (cnt[nb] == 0) killer[nb] = p;
            if (cnt[nb] < 2) ++cnt[nb];
        };

        if (m == 0) {
            for (int nb = 0; nb < 256; ++nb) {
                const bool odd = ((n_odd_ ? 1 : 0) + ((nb - old_b) & 1)) & 1;
                if (!odd) cnt[nb] = 2;
            }
        } else if (!n_odd_) {
            for (int nb = 0; nb < 256; ++nb) cnt[nb] = 2;
        }

        const std::uint32_t hi = static_cast<std::uint32_t>(255 - old_b);
        const std::size_t nb_cnt = pv_.big.size();
        for (std::size_t i = 0; i < nb_cnt; ++i) {
            const std::uint32_t p = pv_.big[i];
            const std::uint32_t t = t_big_[i];
            const std::uint32_t d = t ? p - t : 0;
            const std::uint32_t back = p - d;
            if (d <= hi) {
                mark(static_cast<int>(old_b) + static_cast<int>(d), p);
            } else if (back <= old_b) {
                mark(static_cast<int>(old_b) - static_cast<int>(back), p);
            }
            const std::uint32_t k = static_cast<std::uint8_t>(0u - t * pv_.big_inv[i]);
            t_big_[i] = static_cast<std::uint32_t>((static_cast<std::uint64_t>(t) +
                                                    static_cast<std::uint64_t>(k) * p) >> 8);
        }

        const std::size_t ns_cnt = pv_.small.size();
        for (std::size_t i = 0; i < ns_cnt; ++i) {
            const std::uint32_t p = pv_.small[i];
            const std::uint32_t t = t_small_[i];
            const std::uint32_t d = t ? p - t : 0;
            std::int64_t v = static_cast<std::int64_t>(d) - static_cast<std::int64_t>(p);
            while (v > -255) v -= static_cast<std::int64_t>(p);
            for (; v <= 255; v += static_cast<std::int64_t>(p)) {
                if (v < -255) continue;
                const int nb = static_cast<int>(old_b) + static_cast<int>(v);
                if (nb >= 0 && nb < 256) mark(nb, p);
            }
            const std::uint32_t k = static_cast<std::uint8_t>(0u - t * pv_.small_inv[i]);
            t_small_[i] = static_cast<std::uint32_t>((static_cast<std::uint64_t>(t) +
                                                      static_cast<std::uint64_t>(k) * p) >> 8);
        }
    }

private:
    prime_view pv_;
    std::vector<std::uint32_t> t_big_;
    std::array<std::uint32_t, 64> t_small_{}; 
    bool n_odd_ = true;
};

struct job {
    std::uint32_t pos;
    std::uint8_t nb;
    std::uint32_t p;
};

inline result search(const mpz_class& n, const mpz_class& e) {
    const long long t0 = now_ns();
    const unsigned int nbits = mpz_sizeinbase(n.get_mpz_t(), 2);
    const std::vector<unsigned char> bytes = to_bytes_be(n);
    const std::size_t nbytes = bytes.size();
    if (nbytes == 0) throw std::invalid_argument("one_byte_patchN: n must be positive");

    const bool n_odd = mpz_odd_p(n.get_mpz_t()) != 0;

    std::vector<mpz_class> powers(nbytes);
    powers[0] = 1;
    for (std::size_t m = 1; m < nbytes; ++m)
        mpz_mul_2exp(powers[m].get_mpz_t(), powers[m - 1].get_mpz_t(), 8);

    const std::size_t nthreads = std::max(1u, std::thread::hardware_concurrency());
    const std::uint32_t limit = sieve_limit_for(nbits);

    struct pick {
        bool found = false;
        std::size_t pos = 0;
        std::uint8_t nb = 0;
        mpz_class np, d;
    };
    pick best, cur;
    std::mutex res_mu;

    std::atomic<std::size_t> pos_left{0};

    sieve sv;
    sv.init(n, limit, n_odd);

    std::mutex qmu;
    std::condition_variable qcv;
    std::deque<job> q;
    bool prod_done = false;

    auto producer = [&]() {
        std::uint8_t cnt[256];
        std::uint32_t killer[256];
        std::array<job, 256> batch;
        for (std::size_t m = 0; m < nbytes; ++m) {
            const std::size_t pos = nbytes - 1 - m;
            const std::uint8_t old_b = bytes[pos];
            sv.scan_position(m, old_b, cnt, killer);
            std::size_t n_batch = 0;
            for (int nb = 0; nb < 256; ++nb) {
                if (nb == old_b || (pos == 0 && nb == 0)) continue;
                if (cnt[nb] != 1) continue;
                batch[n_batch++] = job{static_cast<std::uint32_t>(pos),
                                       static_cast<std::uint8_t>(nb), killer[nb]};
            }
            if (n_batch == 0) continue;
            {
                std::lock_guard<std::mutex> lk(res_mu);
                cur.found = false;
            }
            {
                std::lock_guard<std::mutex> lk(qmu);
                pos_left.store(n_batch);
                for (std::size_t i = 0; i < n_batch; ++i) q.push_back(batch[i]);
            }
            qcv.notify_all();
            {
                std::unique_lock<std::mutex> lk(qmu);
                qcv.wait(lk, [&] { return pos_left.load() == 0; });
            }
            {
                std::lock_guard<std::mutex> lk(res_mu);
                if (cur.found) {
                    best = cur;
                    break;
                }
            }
        }
        {
            std::lock_guard<std::mutex> lk(qmu);
            prod_done = true;
        }
        qcv.notify_all();
    };

    auto consumer = [&]() {
        mpz_class np, q_, phi, d, qm1, gg;
        while (true) {
            job j;
            {
                std::unique_lock<std::mutex> lk(qmu);
                qcv.wait(lk, [&] { return !q.empty() || prod_done; });
                if (q.empty()) {
                    if (prod_done) return;
                    continue;
                }
                j = q.front();
                q.pop_front();
            }

            {
                const std::size_t m = nbytes - 1 - j.pos;
                const int delta = static_cast<int>(j.nb) - static_cast<int>(bytes[j.pos]);
                mpz_set(np.get_mpz_t(), n.get_mpz_t());
                if (delta > 0) mpz_addmul_ui(np.get_mpz_t(), powers[m].get_mpz_t(), (unsigned long)delta);
                else mpz_submul_ui(np.get_mpz_t(), powers[m].get_mpz_t(), (unsigned long)-delta);
                mpz_divexact_ui(q_.get_mpz_t(), np.get_mpz_t(), j.p);

                bool coprime = mpz_gcd_ui(nullptr, e.get_mpz_t(),
                                          (unsigned long)(j.p - 1)) == 1;
                if (coprime) {
                    mpz_sub_ui(qm1.get_mpz_t(), q_.get_mpz_t(), 1);
                    if (mpz_fits_ulong_p(e.get_mpz_t())) {
                        coprime = mpz_gcd_ui(nullptr, qm1.get_mpz_t(),
                                             mpz_get_ui(e.get_mpz_t())) == 1;
                    } else {
                        mpz_gcd(gg.get_mpz_t(), e.get_mpz_t(), qm1.get_mpz_t());
                        coprime = mpz_cmp_ui(gg.get_mpz_t(), 1) == 0;
                    }
                }

                if (coprime) {
                    if (mpz_probab_prime_p(q_.get_mpz_t(), 4) != 0 && q_ != j.p) {
                        phi = (mpz_class(j.p) - 1) * (q_ - 1);
                        mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi.get_mpz_t());
                        std::lock_guard<std::mutex> lk(res_mu);
                        if (!cur.found || j.nb < cur.nb) {
                            cur.found = true;
                            cur.pos = j.pos;
                            cur.nb = j.nb;
                            cur.np = np;
                            cur.d = d;
                        }
                    }
                }
            }
            if (pos_left.fetch_sub(1, std::memory_order_acq_rel) == 1) qcv.notify_all();
        }
    };

    std::thread prod(producer);
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (std::size_t i = 0; i < nthreads; ++i) pool.emplace_back(consumer);
    prod.join();
    for (auto& th : pool) th.join();

    result r;
    r.found = best.found;
    r.e = e;
    r.seconds = (now_ns() - t0) / 1e9;
    if (best.found) {
        r.byte_index = best.pos;
        r.new_byte = best.nb;
        r.new_n = best.np;
        r.new_d = best.d;
    }
    return r;
}

}  // namespace detail

inline result find(const mpz_class& n, const mpz_class& e = mpz_class(65537L)) {
    if (n <= 0) throw std::invalid_argument("one_byte_patchN: n must be positive");
    if (e < 3 || mpz_even_p(e.get_mpz_t()))
        throw std::invalid_argument("one_byte_patchN: e must be an odd integer >= 3");
    return detail::search(n, e);
}

}  // namespace one_byte_patchN

#endif  // ONE_BYTE_PATCHN_HPP
