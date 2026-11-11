
# Scalable Bloom Filter

A Scalable Bloom Filter (SBF) built on top of a classic Bloom Filter.
It grows as a linked list of filters with geometrically decreasing
target FPRs, so the *overall* false-positive rate stays bounded by a
single user-supplied $P$ no matter how many elements are inserted.

---

## Table of Contents

- [Features](#features)
- [Background](#background)
- [API Overview](#api-overview)
- [Quick Start](#quick-start)
- [Build](#build)
- [Parameter Allocation](#parameter-allocation)
- [Design Notes](#design-notes)
- [Thread Safety](#thread-safety)
- [Testing](#testing)
- [License](#license)

---

## Features

- **Bounded FPR under growth**: total FPR $\le P$ for any insertion
  count, by allocating block $i$ a target of $P(1-r)r^i$ with
  $r = 0.5$.
- **Doubly-linked blocks**: newest-first query (`cur` -> `head` via
  `prev`) hits the most likely block first, and avoids querying old
  saturated blocks when possible.
- **Automatic expansion**: when the current block is full, a new block
  is appended with capacity $\times 2$ and target FPR $\times 0.5$.
- **Arbitrary byte keys**: keys are `(const void* key, uint64_t len)`,
  not NUL-terminated strings. Binary keys with embedded `'\0'` are
  supported, matching the underlying Bloom Filter contract.
- **Introspection**: block count, element count, current FPR, target
  FPR.
- **`Clear` without reallocation**: resets all bit arrays and counters,
  keeps the block structure, and resets `cur = head` so future inserts
  go back to the first block.
- **Defensive**: every public function handles `NULL`; invalid $P$
  rejected; `MAX_BLOCK_COUNT` caps unbounded chain growth.

---

## Background

A classic Bloom Filter cannot grow: once the bit array fills up, the
FPR climbs past the target. A Scalable Bloom Filter (Almeida et al.,
2007) solves this by chaining multiple filters.

When block $i$ reaches capacity, a new block $i+1$ is created with

$$
C_{i+1} = 2 C_i, \qquad p_{i+1} = r \, p_i \quad (r = 0.5 \text{ by default})
$$

The initial block uses

$$
C_0 = \text{DEFAULT\_EXPECTED\_ELEMENT\_CAPACITY} = 10^6, \qquad
p_0 = P(1 - r) = P/2
$$

Because the target FPRs form a geometric series,

$$
\sum_{i \ge 0} p_0 r^i
  = \frac{p_0}{1 - r}
  = \frac{P(1-r)}{1-r}
  = P
$$

the **union** false-positive rate over all blocks is bounded by $P$.

---

## API Overview

```c
#include "scalable_bloom_filter.h"

/* Memory management */
SBF* SBF_Create(double p);
void SBF_Clear(SBF* sbf);
void SBF_Destroy(SBF* sbf);

/* Operating */
bool SBF_Element_Mark(SBF* sbf, const void* key, uint64_t len);
bool SBF_Element_Query(const SBF* const sbf, const void* key, uint64_t len);

/* Accessors */
uint64_t SBF_Block_Count(const SBF* const sbf);
uint64_t SBF_Ele_Count(const SBF* const sbf);
double   SBF_CFPR(const SBF* const sbf);   /* current overall FPR */
double   SBF_TFPR(const SBF* const sbf);   /* target overall FPR */
```

### Return-value conventions

| Function | Returns | On `NULL` / invalid input |
|---|---|---|
| `Create` | new SBF or `NULL` | `NULL` |
| `Mark` / `Query` | `true` / `false` | `false` |
| `Block_Count` / `Ele_Count` | count | `UINT64_MAX` |
| `CFPR` / `TFPR` | probability in $[0,1]$ | `NAN` |

### Key semantics for `Mark` / `Query`

`key` is treated as a raw byte sequence of `len` bytes. It does **not**
need to be NUL-terminated; embedded `'\0'` bytes are allowed.

- `key != NULL`, `len >= 0`: normal case.
- `key == NULL`, `len == 0`: valid "empty key"; `MurmurHash3_x64_128(NULL, 0, ...)`
  is well-defined, and the filter accepts it. `Mark(sbf, NULL, 0)` succeeds,
  and `Query(sbf, NULL, 0)` returns `true` afterward.
- `key == NULL`, `len > 0`: invalid; both functions return `false`.

This matches the contract of `bloom_filter.h`.

### Notes on `Mark` returning `false`

`SBF_Element_Mark` returns `false` in these cases:

1. `sbf` is `NULL`, or `key == NULL && len > 0`.
2. The current block is full and expansion failed, for any of:
   - the block count has reached `MAX_BLOCK_COUNT`;
   - capacity growth would overflow `uint64_t`;
   - the next target FPR would leave the open interval $(0, 1)$;
   - allocating the new block failed (including `Bloom_Filter_Create`
     failing because the recomputed bit size is too large).

It does **not** mean "element already present". To test membership, use
`SBF_Element_Query`. A failed expansion does not increment
`SBF_Ele_Count`.

---

## Quick Start

```c
#include "scalable_bloom_filter.h"
#include <stdio.h>
#include <string.h>

int main(void) {
    /* overall target FPR 1% */
    SBF* sbf = SBF_Create(0.01);
    if (sbf == NULL) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }

    for (int i = 0; i < 2500000; i++) {
        char key[32];
        int n = snprintf(key, sizeof(key), "user_%d", i);
        SBF_Element_Mark(sbf, key, (uint64_t)n);
    }

    printf("blocks      = %llu\n", (unsigned long long)SBF_Block_Count(sbf));
    printf("elements    = %llu\n", (unsigned long long)SBF_Ele_Count(sbf));
    printf("target FPR  = %.6f\n", SBF_TFPR(sbf));
    printf("current FPR = %.6f\n", SBF_CFPR(sbf));

    SBF_Destroy(sbf);
    return 0;
}
```

Compile (depends on `bloom_filter` and `murmur3`):

```bash
gcc -std=c11 -O2 -Wall -Wextra -Werror \
    bloom_filter.c murmur3.c scalable_bloom_filter.c main.c -lm -o demo
./demo
```

$2.5 \times 10^6$ inserts trigger at least two expansions; `blocks`
will be $\ge 3$.

---

## Build

Depends on `bloom_filter.c` / `bloom_filter.h` and `murmur3.c` /
`murmur3.h`. Compile all together:

```bash
# objects
gcc -std=c11 -O2 -Wall -Wextra -Werror -c bloom_filter.c
gcc -std=c11 -O2 -Wall -Wextra -Werror -c murmur3.c
gcc -std=c11 -O2 -Wall -Wextra -Werror -c scalable_bloom_filter.c

# static library
ar rcs libsbf.a bloom_filter.o murmur3.o scalable_bloom_filter.o
```

Suggested flags: `-std=c11 -O2 -Wall -Wextra -Werror`. Warning-clean
under GCC and Clang.

### Tunable macros

Defined at the top of `scalable_bloom_filter.c`, overridable via `-D`:

| Macro | Default | Meaning |
|---|---|---|
| `DEFAULT_EXPECTED_ELEMENT_CAPACITY` | `1000000` | initial block capacity |
| `CAPACITY_EXPANSION_FACTOR` | `2` | capacity growth per block |
| `FPR_SHRINK_FACTOR` | `0.5` | target-FPR shrink per block |
| `MAX_BLOCK_COUNT` | `64` | hard cap on chain length |

Example:

```bash
gcc -std=c11 -O2 -DDEFAULT_EXPECTED_ELEMENT_CAPACITY=100000 \
    bloom_filter.c murmur3.c scalable_bloom_filter.c main.c -lm -o demo
```

---

## Parameter Allocation

Let $P$ be the user-supplied overall target FPR and $r =
\text{FPR\_SHRINK\_FACTOR}$.

| Block $i$ | Capacity $C_i$ | Target FPR $p_i$ |
|---|---|---|
| $0$ | $C$ | $P(1-r)$ |
| $1$ | $2C$ | $P(1-r)r$ |
| $2$ | $4C$ | $P(1-r)r^2$ |
| $i$ | $C \cdot 2^i$ | $P(1-r)r^i$ |

Sum of target FPRs:

$$
\sum_{i \ge 0} P(1-r) r^i
  = \frac{P(1-r)}{1-r}
  = P
$$

So the **union bound** on false positives is exactly $P$.

The `CFPR` accessor computes the **current** aggregate estimate using
each block's live element count:

$$
\mathrm{CFPR} = 1 - \prod_{i} \left(1 - p_i\right),
\qquad
p_i = \left(1 - e^{-k_i n_i / m_i}\right)^{k_i}
$$

where $n_i$ is the **actual** number of elements currently in block
$i$, not its design capacity. Because $n_i \le C_i$,
$\mathrm{CFPR} \le \mathrm{TFPR}$. The product is evaluated in log
space (`-expm1(sum log1p(-p_i))`) to avoid underflow and cancellation
when many blocks accumulate.

---

## Design Notes

### Doubly-linked blocks

```
head <-> ... <-> cur
        prev  next
```

`next` points toward newer blocks, `prev` toward older blocks.
`SBF_Element_Query` walks from `cur` backward via `prev`:

- newest block is checked first (highest hit probability),
- older, saturated blocks (highest FPR) are checked last.

This is both faster and more accurate than a forward-only walk from
`head`.

### Why not "reverse the list"?

Insertion is always at the tail (`cur`). A doubly-linked list gives
$O(1)$ append **and** $O(1)$ backward query, without having to reverse
or re-index anything.

### `Clear` does not shrink

`SBF_Clear` zeroes all bit arrays, resets counters, and sets
`cur = head`. It does **not** free expanded blocks, does **not** change
`sbf_block_count`, and does **not** change `sbf_tfpr` (both are design
constants / structural state). This avoids repeated `malloc`/`free`
cycles when a workload alternates between fill-and-clear phases.

If you need to return to a single-block state, `Destroy` and call
`Create` again.

### Expansion failure is recoverable, but not rollback-safe

`SBF_Expand` refuses to grow (and `SBF_Element_Mark` returns `false`)
when any of the following holds:

- the chain already has `MAX_BLOCK_COUNT` blocks;
- `prev_capacity * CAPACITY_EXPANSION_FACTOR` would overflow `uint64_t`;
- the next target FPR would leave the open interval $(0, 1)$;
- allocating the new block fails.

In every failure case the chain is left unchanged and the caller can
retry later. However, if an expansion **succeeds** but the subsequent
insert fails, the newly appended (empty) block stays in the chain. This
is harmless for correctness (queries simply find nothing there) but
means `SBF_Block_Count` can grow even when `Mark` returns `false`.

### `CFPR` is a live estimate, `TFPR` is a design bound

`SBF_CFPR` calls `Bloom_Filter_CFPR` on each block, using that block's
**current** `ele_count`. It answers "if I queried a non-member right
now, what is the probability that at least one block returns true?".

`SBF_TFPR` returns the user-supplied $P$, which is the **union bound**
under the assumption that every block is filled to its design
capacity. Since actual $n_i \le C_i$, `CFPR <= TFPR` in practice.

---

## Thread Safety

**Not thread-safe.** Concurrent `Mark` on the same instance races on
the bit array, on the element counter, and on `cur` during expansion.
Use one SBF per thread, or guard shared instances with a mutex.

`SBF_Element_Query` is read-only, but is only safe concurrently with
other readers if no writer is active.

---

## Testing

Two test suites are provided. Both are compiled with
`-std=c11 -O2 -Wall -Wextra -Werror` and run in CI with
AddressSanitizer + UndefinedBehaviorSanitizer.

```bash
# unit tests for the classic Bloom Filter
gcc -std=c11 -O2 -Wall -Wextra -Werror \
    bloom_filter.c murmur3.c bf_test.c -lm -o test_bf
./test_bf

# unit + integration tests for the scalable filter
gcc -std=c11 -O2 -Wall -Wextra -Werror \
    bloom_filter.c murmur3.c scalable_bloom_filter.c sbf_test.c -lm -o test_sbf
./test_sbf
```

Coverage for the SBF suite:

- Constructor validation ($P \le 0$, $P \ge 1$)
- `NULL` defense for every public function, including the
  `key == NULL && len == 0` (valid empty key) vs
  `key == NULL && len > 0` (invalid) distinction
- $\mathrm{CFPR} = 0$ on an empty filter; $\mathrm{TFPR}$ matches the
  input
- **No false negatives across expansions** ($1.5 \times 10^6$ inserts
  -> 2 blocks), with `Ele_Count` matching the insert count
- **`Clear` resets state** across multiple blocks: block count is
  preserved, `TFPR` is unchanged, all previously inserted keys miss,
  and the next insert goes into `head`
- **Realistic scenario**: 800k "registered users", 400k "unknown
  users"; no false negatives, measured FP rate within `0.005` absolute
  distance of `CFPR`
- **Expansion stress**: $3.5 \times 10^6$ inserts -> at least 3 blocks,
  sampled queries across all keys all succeed (exercises the `prev`
  chain)

---

## License

MIT.