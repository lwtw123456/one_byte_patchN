// demo: thin CLI wrapper around the library (the library itself does no I/O).
//   demo <hex modulus> [e]
// All output is ASCII, so it renders correctly in a GBK console.
#include "one_byte_patchN.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static std::string hex_up(const mpz_class& v) {
    char* s = mpz_get_str(nullptr, 16, v.get_mpz_t());
    std::string r(s);
    void (*freefn)(void*, std::size_t) = nullptr;
    mp_get_memory_functions(nullptr, nullptr, &freefn);
    freefn(s, r.size() + 1);
    for (char& c : r) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return r;
}

static std::string dec(const mpz_class& v) {
    char* s = mpz_get_str(nullptr, 10, v.get_mpz_t());
    std::string r(s);
    void (*freefn)(void*, std::size_t) = nullptr;
    mp_get_memory_functions(nullptr, nullptr, &freefn);
    freefn(s, r.size() + 1);
    return r;
}

static std::vector<unsigned char> bytes_be(const mpz_class& v) {
    std::vector<unsigned char> out;
    out.resize((mpz_sizeinbase(v.get_mpz_t(), 2) + 7) / 8);
    std::size_t count = 0;
    mpz_export(out.data(), &count, 1, 1, 1, 0, v.get_mpz_t());
    out.resize(count);
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::fprintf(stderr, "usage: %s <hex modulus> [e]\n", argv[0]);
        return 2;
    }
    mpz_class n(argv[1], 16);
    mpz_class e(65537L);
    if (argc == 3) mpz_set_str(e.get_mpz_t(), argv[2], 0);

    try {
        const auto r = one_byte_patchN::find(n, e);
        if (!r.found) {
            std::printf("not found (n = %lu bit, e = %s, %.6f s)\n",
                        (unsigned long)mpz_sizeinbase(n.get_mpz_t(), 2), dec(r.e).c_str(), r.seconds);
            return 1;
        }
        const std::vector<unsigned char> nb = bytes_be(n);
        const unsigned char old_byte = r.byte_index < nb.size() ? nb[r.byte_index] : 0;
        std::printf("e: %s\n", dec(r.e).c_str());
        std::printf("byte: index %zu (depth %zu), %02X -> %02X\n", r.byte_index,
                    nb.size() - 1 - r.byte_index, old_byte, r.new_byte);
        std::printf("new_n: %s\n", hex_up(r.new_n).c_str());
        std::printf("new_d: %s\n", hex_up(r.new_d).c_str());
        std::printf("time: %.6f s\n", r.seconds);
        return 0;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "error: %s\n", ex.what());
        return 2;
    }
}
