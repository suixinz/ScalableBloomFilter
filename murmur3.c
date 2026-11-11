
#include "murmur3.h"
#include <string.h> //memcpy()


/*
 * rotl64 - rotate a 64-bit value left by r bits.
 *
 * CONTRACT (caller must guarantee): 1 <= r <= 63.
 * No runtime check is performed, by design: every call site in this
 * file passes a compile-time constant in [1, 63], so the contract is
 * satisfied by construction and a branch would be pure overhead.
 *
 * If you modify this function or add a new call site, be aware that
 * r outside [1, 63] is undefined behavior, because the body performs
 * two shifts with different UB conditions (C11 6.5.7p3: the shift
 * count must be >= 0 and < 64):
 *
 *   - x << r         requires 0 <= r < 64
 *   - x >> (64 - r)  requires 0 <= (64 - r) < 64
 *                    i.e.    1 <= r <= 64
 *
 * Intersecting the two gives 1 <= r <= 63. In particular:
 *
 *   - r == 0   : x << 0 is fine, but 64 - r == 64, so x >> 64 has a
 *                shift count equal to the type width -> UB.
 *   - r == 64  : x << 64 has a shift count equal to the type width
 *                -> UB; x >> 0 is fine, but the left shift already
 *                makes the whole expression UB.
 *   - r < 0    : both shifts have a negative count -> UB.
 *   - r > 64   : x << r has count >= 64, and 64 - r is negative
 *                -> UB on both sides.
 *
 * If you ever need r == 0, or a value that may fall outside [1, 63]
 * (e.g. a rotation amount coming from user input), do NOT relax the
 * contract silently. Replace the body with the masked form, which is
 * well-defined for any unsigned r and matches the mathematical
 * periodicity of rotation (period 64):
 *
 *     static inline uint64_t rotl64(uint64_t x, unsigned r) {
 *         r &= 63;                              // 0 <= r <= 63
 *         return (x << r) | (x >> ((64 - r) & 63));
 *     }
 *
 * Note that the masked form changes the signature (int8_t -> unsigned)
 * and turns every input into well-defined behavior; it also makes
 * r == 64 behave as r == 0, which is the mathematically correct
 * rotation by 64 (identity), not an error.
 */
static inline uint64_t rotl64(uint64_t x, int8_t r){

    return (x << r) | (x >> (64 - r));

}

static inline uint64_t fmix64(uint64_t k){

    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;

    return k;
}

static inline uint64_t getblock64(const uint8_t* p){

    uint64_t v;
    
    memcpy(&v, p, 8);

    return v;

}

void MurmurHash3_x64_128(const void* key, size_t len, uint32_t seed, uint64_t out[2]){

    const uint8_t* data = (const uint8_t*)key;
    const size_t nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;

    /* body */
    const uint8_t* blocks = data;

    for(size_t i = 0; i < nblocks; i++){

        uint64_t k1 = getblock64(blocks + i * 16 + 0);
        uint64_t k2 = getblock64(blocks + i * 16 + 8);

        k1 *= c1; 
        k1 = rotl64(k1, 31); 
        k1 *= c2; 
        h1 ^= k1;
        h1 = rotl64(h1, 27); 
        h1 += h2; 
        h1 = h1 * 5 + 0x52dce729ULL;

        k2 *= c2; 
        k2 = rotl64(k2, 33); 
        k2 *= c1; 
        h2 ^= k2;
        h2 = rotl64(h2, 31); 
        h2 += h1; 
        h2 = h2 * 5 + 0x38495ab5ULL;

    }

    /* tail */
    const uint8_t* tail = data + nblocks * 16;
    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch (len & 15){

        case 15: k2 ^= ((uint64_t)tail[14]) << 48; /* fallthrough */
        case 14: k2 ^= ((uint64_t)tail[13]) << 40; /* fallthrough */
        case 13: k2 ^= ((uint64_t)tail[12]) << 32; /* fallthrough */
        case 12: k2 ^= ((uint64_t)tail[11]) << 24; /* fallthrough */
        case 11: k2 ^= ((uint64_t)tail[10]) << 16; /* fallthrough */
        case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;  /* fallthrough */
        case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;  

                 k2 *= c2; 
                 k2 = rotl64(k2, 33); 
                 k2 *= c1; h2 ^= k2;
                 /* fallthrough */
                 
        case  8: k1 ^= ((uint64_t)tail[ 7]) << 56; /* fallthrough */
        case  7: k1 ^= ((uint64_t)tail[ 6]) << 48; /* fallthrough */
        case  6: k1 ^= ((uint64_t)tail[ 5]) << 40; /* fallthrough */
        case  5: k1 ^= ((uint64_t)tail[ 4]) << 32; /* fallthrough */
        case  4: k1 ^= ((uint64_t)tail[ 3]) << 24; /* fallthrough */
        case  3: k1 ^= ((uint64_t)tail[ 2]) << 16; /* fallthrough */
        case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;  /* fallthrough */
        case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;  

                 k1 *= c1; 
                 k1 = rotl64(k1, 31); 
                 k1 *= c2; h1 ^= k1;
                
        default: break;

    }

    /* finalization */
    h1 ^= (uint64_t)len;
    h2 ^= (uint64_t)len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    out[0] = h1;
    out[1] = h2;

}
