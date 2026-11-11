
# Bloom Filter

A classic Bloom Filter implemented in C11, with a clean O(1) API,
MurmurHash3-based double hashing with coprime step normalization, and
full parameter introspection (target FPR vs. current FPR).

---

## Table of Contents

- [Features](#features)
- [Background](#background)
- [API Overview](#api-overview)
- [Quick Start](#quick-start)
- [Build](#build)
- [Parameter Formulas](#parameter-formulas)
- [Design Notes](#design-notes)
- [Thread Safety](#thread-safety)
- [Testing](#testing)
- [License](#license)

---

## Features

- **Correct math**: $m$ and $k$ derived from $n$ and $p$ using the
  standard formulas; the bit array is aligned up to a byte boundary.
- **Low false-positive clustering**: double hashing ($h_1 + i \cdot h_2$)
  with $h_2$ normalized to be coprime with the bit size via GCD, so the
  $k$ probes never collapse onto a small subset of slots.
- **MurmurHash3 x64 128-bit**: a single 128-bit hash call produces two
  independent 64-bit values, used directly as $h_1$ and $h_2$. No
  separate hash functions, no string-specific assumptions.
- **Arbitrary byte keys**: keys are `(const void* key, uint64_t len)`,
  not NUL-terminated strings. Binary keys with embedded `'\0'` are
  fully supported; `"a\0b"` and `"a\0c"` are distinct keys.
- **Introspection**: query capacity, element count, target FPR, and the
  *current* FPR at any time.
- **Defensive**: every public function accepts `NULL` safely; invalid
  constructor parameters return `NULL`; overflow of $m$ is guarded;
  $n > \mathrm{UINT64\_MAX}/2$ is rejected to keep the element counter
  from wrapping.
- **No global state**: multiple independent filters can coexist.

---

## Background

A Bloom Filter answers "is this element possibly in the set?" in
$O(k)$ time and $O(m)$ bits, with a tunable false-positive rate $p$
and no false negatives.

Given an expected element count $n$ and a target FPR $p \in (0,1)$:

$$
m = \left\lceil \frac{-n \ln p}{(\ln 2)^2} \right\rceil
$$

$$
k = \max\!\left( \left\lceil \frac{-\ln p}{\ln 2} \right\rceil,\ 1 \right)
$$

The API is intentionally minimal: **create → mark → query → destroy**,
plus read-only accessors.

---

## API Overview

```c
#include "bloom_filter.h"

/* Memory management */
Bloom_Filter* Bloom_Filter_Create(uint64_t n, double p);
void          Bloom_Filter_Clear(Bloom_Filter* bf);
void          Bloom_Filter_Destroy(Bloom_Filter* bf);

/* Operating */
bool Bloom_Filter_Element_Mark(Bloom_Filter* bf, const void* key, uint64_t len);
bool Bloom_Filter_Element_Query(const Bloom_Filter* const bf, const void* key, uint64_t len);

/* Accessors */
uint64_t Bloom_Filter_Element_Count(const Bloom_Filter* const bf);
uint64_t Bloom_Filter_Capacity(const Bloom_Filter* const bf);
double   Bloom_Filter_TFPR(const Bloom_Filter* const bf);   /* target FPR */
double   Bloom_Filter_CFPR(const Bloom_Filter* const bf);   /* current FPR */
bool     Bloom_Filter_IsFull(const Bloom_Filter* const bf);
```

### Return-value conventions

| Function | Returns | On `NULL` / invalid input |
|---|---|---|
| `Create` | new filter or `NULL` | `NULL` |
| `Mark` / `Query` | `true` / `false` | `false` |
| `Element_Count` / `Capacity` | count / capacity | `UINT64_MAX` |
| `TFPR` / `CFPR` | probability in $[0,1]$ | `NAN` |
| `IsFull` | full flag | `true` |

### Key semantics for `Mark` / `Query`

`key` is treated as a raw byte sequence of `len` bytes. It does **not**
need to be NUL-terminated; embedded `'\0'` bytes are allowed.

- `key != NULL`, `len >= 0`: normal case.
- `key == NULL`, `len == 0`: valid "empty key"; `MurmurHash3_x64_128(NULL, 0, ...)`
  is well-defined, and the filter accepts it. `Mark(bf, NULL, 0)` succeeds,
  and `Query(bf, NULL, 0)` returns `true` afterward.
- `key == NULL`, `len > 0`: invalid; both functions return `false`.

---

## Quick Start

```c
#include "bloom_filter.h"
#include <stdio.h>
#include <string.h>

int main(void) {
    /* capacity 1,000,000, target FPR 1% */
    Bloom_Filter* bf = Bloom_Filter_Create(1000000, 0.01);
    if (bf == NULL) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }

    const char* alice = "alice";
    const char* bob   = "bob";

    Bloom_Filter_Element_Mark(bf, alice, strlen(alice));
    Bloom_Filter_Element_Mark(bf, bob,   strlen(bob));

    printf("alice: %s\n", Bloom_Filter_Element_Query(bf, alice, strlen(alice)) ? "maybe yes" : "no");
    printf("carol: %s\n", Bloom_Filter_Element_Query(bf, "carol", 5) ? "maybe yes" : "no");

    printf("target FPR  = %.6f\n", Bloom_Filter_TFPR(bf));
    printf("current FPR = %.6f\n", Bloom_Filter_CFPR(bf));

    Bloom_Filter_Destroy(bf);
    return 0;
}
```

Compile:

```bash
gcc -std=c11 -O2 -Wall -Wextra -Werror \
    bloom_filter.c murmur3.c main.c -lm -o demo
./demo
```

---

## Build

The library is two `.c` files (`bloom_filter.c`, `murmur3.c`) plus two
headers (`bloom_filter.h`, `murmur3.h`). No external dependencies
beyond the C standard library and `libm`.

```bash
# object files
gcc -std=c11 -O2 -Wall -Wextra -Werror -c bloom_filter.c murmur3.c

# static library (optional)
ar rcs libbloom_filter.a bloom_filter.o murmur3.o
```

Suggested flags: `-std=c11 -O2 -Wall -Wextra -Werror`. The code is
warning-clean under GCC and Clang.

---

## Parameter Formulas

| Symbol | Meaning | Formula |
|---|---|---|
| $n$ | expected elements | user input, $0 < n \le \mathrm{UINT64\_MAX}/2$ |
| $p$ | target FPR | user input, $0 < p < 1$ |
| $m$ | bit count | $m = \left\lceil -n \ln p / (\ln 2)^2 \right\rceil$ |
| $k$ | hash functions | $k = \max\!\left(\left\lceil -\ln p / \ln 2 \right\rceil,\ 1\right)$ |
| $B$ | allocated bits | $B = \mathrm{AlignUp}(m, 8)$ |
| $\mathrm{CFPR}$ | current FPR | $\mathrm{CFPR} = \left(1 - e^{-k n_{\mathrm{ins}} / B}\right)^{k}$ |

Example: $n = 10^6$, $p = 0.01$ gives $m \approx 9.59 \times 10^6$
bits (~1.2 MB) and $k = 7$.

---

## Design Notes

### MurmurHash3 x64 128-bit

The filter uses `MurmurHash3_x64_128(key, len, seed=0, out[2])` to
produce $h_1 = \text{out}[0]$ and $h_2 = \text{out}[1]$. The two
64-bit halves are statistically independent, so a single hash call
suffices for the Kirsch–Mitzenmacher double-hashing scheme.

Earlier versions used BKDR + DJB with a Murmur3 finalizer; that worked
but the two hashes were correlated, and the FPR silently degraded on
adversarial or structured keys. Murmur3 is the standard choice: fast,
well-distributed, and has published test vectors that `bf_test.c`
verifies against.

### Double hashing with coprime step

Instead of computing $k$ independent hashes per element, the filter
uses

$$
\mathrm{pos}_i = (h_1 + i \cdot h_2) \bmod B, \qquad i = 0, 1, \dots, k-1
$$

where $B$ is the bit size. This is the Kirsch–Mitzenmacher
double-hashing scheme. To prevent the probes from collapsing onto a few
slots (which happens when $\gcd(h_2, B) > 1$), $h_2$ is normalized:

1. Clamp to $[1, B-1]$ — a step of $0$ would repeat the same bit.
2. Increment until $\gcd(h_2, B) = 1$.

This guarantees the $k$ probes traverse $k$ distinct residues modulo
$B$, matching the FPR bound of $k$ truly independent hashes.

### Why cast before negation

```c
double m_d = ceil(-(double)n * ln_p / (ln_2 * ln_2));   /* correct */
/* double m_d = ceil(-n * ln_p / ...);                   WRONG on unsigned n */
```

`n` is `uint64_t`. Applying unary `-` to an unsigned integer wraps
around, so $-n$ becomes $\mathrm{UINT64\_MAX} - n + 1$ — a huge
positive number. Always cast to `double` first.

### Why `n > UINT64_MAX/2` is rejected

`ele_count` is a `uint64_t` and is incremented once per `Mark` call,
including duplicate Marks. The worst case is $2n$ (every element
marked twice). Keeping $n \le \mathrm{UINT64\_MAX}/2$ guarantees
`ele_count` never wraps, which in turn keeps `IsFull` and `CFPR`
reliable.

In practice $n$ is far below this bound; the check exists to make the
contract explicit, not because the case is reachable.

### `IsFull` does not block `Mark`

The classic filter is a pure data structure: `Mark` always writes.
`IsFull` reports $n_{\mathrm{ins}} \ge \text{capacity}$ so callers can
decide when to stop. The scalable filter builds on this convention.

---

## Thread Safety

**Not thread-safe.** Concurrent `Mark` on the same instance races on
`bit_array |= ...` and on the element counter. Use one filter per
thread, or guard shared instances with a mutex.

---

## Testing

A unit-test suite is provided in `bf_test.c`.

```bash
gcc -std=c11 -O2 -Wall -Wextra -Werror \
    bloom_filter.c murmur3.c bf_test.c -lm -o test_bf
./test_bf
```

Coverage:

- **MurmurHash3 known vectors**: `hash("", 0, seed=0) == {0, 0}`,
  plus determinism, seed sensitivity, and independence of the two
  64-bit halves. This is the test that catches a broken hash
  implementation, which would otherwise pass all other tests while
  silently degrading the FPR.
- Constructor validation ($n = 0$, $n > \mathrm{UINT64\_MAX}/2$,
  $p \le 0$, $p \ge 1$)
- `NULL` defense for every public function, including the
  `key == NULL && len == 0` / `key == NULL && len > 0` distinction
- Tiny capacities ($n = 1, 2, 3, 4, 8$), exercising `Normalize_Step`
  on small `bit_size` values
- Binary keys with embedded `'\0'`, including keys that differ only in
  a byte after the first `'\0'`
- No false negatives over 200k inserts
- $\mathrm{CFPR} \approx \mathrm{TFPR}$ when filled to capacity
  (30% relative tolerance)
- Measured FP rate within an absolute distance of $\mathrm{CFPR}$
  (0.005), which is stricter than "within 3× of target" and catches
  hash-quality regressions
- Idempotent `Mark`, and `Clear` resetting all state
- `IsFull` semantics and accessor values

---

## License

MIT.
