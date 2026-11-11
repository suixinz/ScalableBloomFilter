
#include "bloom_filter.h"

#include <stdlib.h>     /**< calloc() / free() */
#include <math.h>       //log() ceil() pow() exp()
#include <string.h>     //memset()
#include "murmur3.h"    //MurmurHash3_x64_128()

#ifndef BITS_PER_BYTE

#define BITS_PER_BYTE ((uint64_t)8)

#endif

struct Bloom_Filter{

    uint64_t        capacity;
    double          tfpr;

    unsigned char*  bit_array;
    uint64_t        bit_size;
    uint64_t        hash_count;

    uint64_t        ele_count;

};

/*

n-->​	Expected number of elements​ to be inserted into the filter
p-->​	Target false positive rate​ (acceptable probability of false positives, between 0 and 1)
m​-->	Size of the bit array​ (total number of bits)
k​-->	Number of hash functions​ used per element

relation: n, p ==> m, k

m = - n * ln(p) / (ln(2))^2;
k = - ln(p) / ln(2);

*/

//==================== Memory Management ====================

static inline uint64_t Align_Up(uint64_t n){
    return (n + BITS_PER_BYTE - 1) & (~(BITS_PER_BYTE - 1));
}

static inline uint64_t Max_U64(uint64_t a, uint64_t b){
    return a > b ? a : b;
}

Bloom_Filter* Bloom_Filter_Create(uint64_t n, double p){

    Bloom_Filter* bf = NULL;

    /*
     * Reject n that could cause ele_count to wrap. ele_count is
     * incremented once per Mark call, including duplicate Marks, so
     * the worst case is 2 * n (if every element is marked twice).
     * Keeping n <= UINT64_MAX / 2 guarantees ele_count never wraps.
     *
     * In practice n is far below this bound; the check exists to make
     * the contract explicit, not because the case is reachable.
     */

    if(n == 0 || n > (UINT64_MAX / 2) || !(p > 0.0 && p < 1.0)){
        goto cleanup;
    }

    bf = (Bloom_Filter*)calloc(1, sizeof(Bloom_Filter));

    if(bf == NULL){
        goto cleanup;
    }

    const double ln_2 = log(2);
    const double ln_p = log(p);

    double m_d = ceil(-(double)n * ln_p / (ln_2 * ln_2));

    /*
    * WARNING: Applying unary minus to an unsigned integer (uint64_t here)
    * does NOT produce a negative number due to unsigned wrap-around.
    *
    * For example:
    *   size_t n = 10000;
    *   -n  →  SIZE_MAX - 10000 + 1  (a huge positive value on 64-bit)
    *
    * Always cast to double first before negation:
    *   double result = -(double)n;   // Correct: -10000.0
    * 2026.10.26
    */

    if (m_d > (double)(UINT64_MAX / BITS_PER_BYTE)){
        goto cleanup;
    }

    const uint64_t m = (uint64_t)m_d;
    const uint64_t k = Max_U64((uint64_t)ceil(- ln_p / ln_2), 1);

    bf->bit_size = Align_Up(m);
    bf->bit_array = (unsigned char*)calloc(bf->bit_size / BITS_PER_BYTE, sizeof(unsigned char));

    if(bf->bit_array == NULL){
        goto cleanup;
    }

    bf->hash_count = k;
    bf->ele_count = 0;

    bf->capacity = n;
    bf->tfpr = p;

    return bf;


    cleanup:

        Bloom_Filter_Destroy(bf);
        return NULL;

}

void Bloom_Filter_Clear(Bloom_Filter* bf){

    if(bf == NULL){
        return;
    }

    memset(bf->bit_array, 0, bf->bit_size / BITS_PER_BYTE);

    bf->ele_count = 0;

}

void Bloom_Filter_Destroy(Bloom_Filter* bf){

    if(bf == NULL){
        return;
    }

    free(bf->bit_array);
    free(bf);

}


//==================== Operating Functions ====================

/* Greatest common divisor (Euclidean algorithm) */
static inline uint64_t GCD_U64(uint64_t a, uint64_t b){

    while(b){
        uint64_t t = a % b;
        a = b;
        b = t;
    }

    return a;

}

/*
 * Normalize_Step - turn an arbitrary hash value into a valid step for
 * double hashing: h1 + i * step (mod bit_size).
 *
 * Precondition: bit_size > 0. Bloom_Filter_Create guarantees this by
 * construction (bit_size = Align_Up(m) >= 8 whenever m > 0), so no
 * runtime check is performed here.
 *
 * Postcondition: 1 <= return value < bit_size, and
 * gcd(return value, bit_size) == 1. The gcd condition ensures the k
 * probes h1 + i*step (i = 0..k-1) visit k distinct residues modulo
 * bit_size, matching the false-positive bound of k truly independent
 * hashes.
 *
 * The while loop terminates: step == 1 always satisfies
 * gcd(1, bit_size) == 1, so the wrap-around branch is unreachable in
 * practice. It is kept as a defensive no-op in case the loop body is
 * ever changed.
 */
static inline uint64_t Normalize_Step(uint64_t h, uint64_t bit_size){
    
    uint64_t step = h % bit_size;

    if(step == 0){
        step = 1;
    }

    while(GCD_U64(step, bit_size) != 1){

        step++;
        if(step >= bit_size){
            step = 1;   
            // wrap-around; actually unreachable: step == 1 always satisfies GCD(1, bit_size) == 1, so the loop exits before wrapping.
        }

    }

    return step;

}

bool Bloom_Filter_Element_Mark(Bloom_Filter* bf, const void* key, uint64_t len){

    if(bf == NULL || (key == NULL && len > 0)){
        return false;
    }

    uint64_t out[2];

    MurmurHash3_x64_128(key, len, 0, out);

    const uint64_t h1 = out[0];
    const uint64_t h2 = Normalize_Step(out[1], bf->bit_size);

    for(uint64_t i = 0; i < bf->hash_count; i++){

        uint64_t pos_i = (h1 + i * h2) % bf->bit_size;

        bf->bit_array[pos_i / BITS_PER_BYTE] |= ((uint64_t)1 << (pos_i % BITS_PER_BYTE));

    }

    bf->ele_count++;

    return true;

}

bool Bloom_Filter_Element_Query(const Bloom_Filter* const bf, const void* key, uint64_t len){

    if(bf == NULL || (key == NULL && len > 0)){
        return false;
    }

    uint64_t out[2];
    MurmurHash3_x64_128(key, len, 0, out);

    const uint64_t h1 = out[0];
    const uint64_t h2 = Normalize_Step(out[1], bf->bit_size);

    for(uint64_t i = 0; i < bf->hash_count; i++){

        uint64_t pos_i = (h1 + i * h2) % bf->bit_size;

        if(bf->bit_array[pos_i / BITS_PER_BYTE] & ((uint64_t)1 << (pos_i % BITS_PER_BYTE))){
            continue;
        }

        return false;

    }

    return true;

}


//==================== Parameter Accessor Functions ====================

uint64_t Bloom_Filter_Element_Count(const Bloom_Filter* const bf){

    if(bf == NULL){
        return UINT64_MAX;
    }

    return bf->ele_count;

}

uint64_t Bloom_Filter_Capacity(const Bloom_Filter* const bf){

    if(bf == NULL){
        return UINT64_MAX;
    }

    return bf->capacity;

}

double Bloom_Filter_TFPR(const Bloom_Filter* const bf){

    if(bf == NULL){
        return NAN;
    }

    return bf->tfpr;

}

double Bloom_Filter_CFPR(const Bloom_Filter* const bf){

    if(bf == NULL){
        return NAN;
    }

    // Formula: p = (1 - e^(-kn/m))^k
    // k: number of hash functions
    // n: number of inserted elements
    // m: size of bit array

    double k = bf->hash_count;
    double n = bf->ele_count;
    double m = bf->bit_size;
    
    double exponent = -k * n / m;
    double cur_p = pow(-expm1(exponent), k);
    
    return cur_p;

}

bool Bloom_Filter_IsFull(const Bloom_Filter* const bf){

    if(bf == NULL){
        return true;
    }

    return bf->ele_count >= bf->capacity;

}

