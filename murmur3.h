#pragma once

#include <stdint.h>
#include <stddef.h>


/*
 * MurmurHash3_x64_128 - 128-bit MurmurHash3 (x64 variant).
 *
 * Deterministic: the same (key, len, seed) always produces the same
 * (out[0], out[1]). This is a hash, not a random number generator.
 *
 * Precondition: `key` may be NULL only when `len == 0`.
 *               `out` must point to a writable uint64_t[2].
 *
 * The two outputs are statistically independent enough to be used
 * directly as h1 and h2 for double hashing.
 *
 * If you change the implementation, verify it against the published
 * MurmurHash3 test vectors (e.g. hash("", 0, seed=0) must produce
 * out[0] == 0 and out[1] == 0). A silent regression here degrades
 * the Bloom Filter's false-positive rate without any test failure
 * elsewhere.
 */
void MurmurHash3_x64_128(const void* key, size_t len, uint32_t seed, uint64_t out[2]);
