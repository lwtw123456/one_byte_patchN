# one_byte_patchN.hpp

RSA 加密 1 字节 patchN 的 C++ header-only 库。

## 依赖

| | |
| --- | --- |
| 编译器 | 支持 `__int128` 的 C++20 编译器（GCC / Clang） |
| 库 | GMP 6.x（含 gmpxx） |

**不要用 MSVC。** MSVC 没有 GCC 风格的 `__int128`，GMP 在 x86_64 上拿不到 64×64→128 位的直接中间结果，`mpn` 的乘法和平方会退化成 32 位拼接的 C 代码。

GMP 的 x86_64 汇编路径（`libgmp.a` 里 680 条 `mulq`/`mulx`）才是 Miller-Rabin 那 99.9% 耗时的来源。

## 编译

Windows（MSYS2 UCRT64）：

```bat
g++ -std=c++20 -O3 -DNDEBUG -Iinclude examples\demo.cpp -o demo.exe -static -s -ffunction-sections -fdata-sections -Wl,--gc-sections -lgmpxx -lgmp
```

Linux / macOS：

```bash
g++ -std=c++20 -O3 -DNDEBUG -Iinclude examples/demo.cpp -o demo -s -ffunction-sections -fdata-sections -Wl,--gc-sections -lgmpxx -lgmp
```

装 GMP：

| 平台 | 命令 |
| --- | --- |
| MSYS2 | `pacman -S mingw-w64-ucrt-x86_64-gmp` |
| Debian / Ubuntu | `apt install libgmp-dev` |
| macOS | `brew install gmp` |

## 用法

见 demo 示例。

`e` 是 `find(n, e)` 的第二个参数，默认 `65537`：

```cpp
auto r = one_byte_patchN::find(n);       // e = 65537
auto r = one_byte_patchN::find(n, 3);    // 指定 e
```

`e` 必须是**奇数且 ≥ 3**，否则抛 `std::invalid_argument`。

`e` 越小，预筛越强，但**解也会跟着变**——同一个 `n` 换个 `e` 可能落到更深的字节位置。

这批用例里 1024 bit 上 `e = 3` 更快；2048 bit 以上 depth 0 没有满足 `gcd(3, φ) = 1` 的解，只能往后扫，反而更慢。

「无解」时 `find()` 返回 `found == false`，不抛异常。

demo 里 `e` 是可选的第二个参数，按 `base 0` 解析（十进制 / `0x` 十六进制 / 前导 `0` 八进制）：

```bat
demo.exe <十六进制模数>           :: e = 65537
demo.exe <十六进制模数> 3         :: e = 3
demo.exe <十六进制模数> 0x10001   :: e = 65537
```

## 实测

本机配置：

| | |
| --- | --- |
| CPU | AMD Ryzen 5 5600（6 核 12 线程，3.5 GHz） |
| 内存 | 16 GB |
| 系统 | Windows 11 Pro（build 26200） |
| 编译器 | GCC 16.2.0（MSYS2 UCRT64，`x86_64-w64-mingw32`） |
| GMP | 6.3.0，静态 `libgmp.a`（fat 构建，按 CPUID 选 `mpn` 变体） |
| 编译选项 | `-O3 -DNDEBUG -static -s --gc-sections` |

`tests/data.txt` 里每个尺寸一个用例，
各跑 15 次取中位数（12 线程）：

| 模数 | 耗时 | 改的字节 |
| --- | --- | --- |
| 1024 bit | 4.5 ms | byte 127，8F → 39 |
| 2048 bit | 14.7 ms | byte 255，77 → 75 |
| 3072 bit | 43.6 ms | byte 383，DF → DD |
| 4096 bit | 84.4 ms | byte 511，BD → BF |

时间几乎全在 Miller-Rabin，和「解落在第几个字节」强相关。

### 与其他实现的速度对比

同一台机器，同样的输入，对比闭源实现 rpatch(https://github.com/nszy007/rpatch) ：

| 模数 | 本库 | rpatch.exe |
| --- | --- | --- |
| 1024 bit | 4.5 ms | 19.2 ms |
| 2048 bit | 14.7 ms | 823 ms |
| 3072 bit | 43.6 ms | 4.85 s |
| 4096 bit | 84.4 ms | 26.4 s |

至于飘云阁的 RSA 一字节替换工具 v1.3 以及开源实现 RSA-PATCH-TOOLKIT(https://github.com/FTSTT/RSA-PATCH-TOOLKIT) ，速度差距更大，4096 bit 上千倍起步。

## 为什么快

对固定的字节位置 `m`，`n'(d) = n + d·256^m`（`d` = 新字节 − 旧字节）。对任一奇素数 `r`：

```
n'(d) ≡ 0 (mod r)   ⟺   d ≡ −n · 256^(−m) (mod r)
```

只要知道 `n·256^(−m) mod r`，就能一次反推出这个位置上 256 种取值里哪些被 `r` 整除。所以扫一遍素数表，`cnt[b]`（改成 `b` 后被几个小素数整除，饱和到 2）和 `killer[b]`（其中最小的那个）就全有了，不用对每个候选挨个试除。

- **换位置只花两次 32 位乘法**：`256^(−(m+1)) = 256^(−m) · 256^(−1) (mod r)`。代码取 `k = −t·r^(−1) mod 256`，用 `t ← (t + k·r) / 256` 更新余数——每个（位置, 素数）就是两次 32 位乘法加移位，不用对大整数重新取模。
- **只有 `cnt == 1` 的候选值得判素**：此时 `n'` 恰好有一个不超过上限的小素因子 `p`，商 `q` 与所有小素数互素，判素命中率高；`cnt ≥ 2` 说明被两个以上小素数整除，直接是合数。
- **`gcd(e, φ) ≠ 1` 再筛一道**：算不出 `d` 的候选不进 Miller-Rabin，`e = 3` 时能扔掉约 3/4。
- **剩下的交给 GMP 和 多线程**：256 个候选里通常只剩几十个（约 1/4）真正要做 Miller-Rabin，走 GMP 的 x86_64 汇编（`mulq`/`mulx`）；producer 按字节位置从低到高批量投递，消费者并行做除法、gcd、判素。
- **素数表预计算**：素数（3e4 / 3e5）和每个素数的 `r⁻¹ mod 256` 提前算好，`p > 255` 与 `p ≤ 255` 分段，热循环里没有 size 分支。

净效果：99.9% 的时间在 Miller-Rabin，但需要判素的候选已经被压到 256 个里的 1/4 左右。

## 素数表：编译期还是运行期

默认**运行期**建一次（堆上）：exe 小 206 KB，代价是进程内第一次 `find()` 多 0.05 ms（≤2048 bit）/ 0.58 ms（>2048 bit）。两种方式生成的表逐字节一致。

要编译期建表（表进 `.rdata`，exe 大 206 KB、零建表成本）：

```bash
g++ ... -DBYTE_PATCH_TABLE_AT_COMPILE_TIME=1 ...
```

> Clang 用编译期建表要加 `-fconstexpr-steps=100000000`。

## 目录

```
include/one_byte_patchN.hpp   库本体（单头文件）
examples/demo.cpp           命令行示例
tests/data.txt              样例模数（1024 / 2048 / 3072 / 4096 bit）
```

## 开源协议

本项目基于 MIT License
开源，可自由用于个人或商业用途，但需保留版权声明。
