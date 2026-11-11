
#include "scalable_bloom_filter.h"

#include <stdint.h>     /* uint64_t */
#include <math.h>       /* pow() */
#include <stdlib.h>     /* calloc() / free() */

/*
 * scalable_bloom_filter.c
 *
 * A Scalable Bloom Filter is a linked list of classic Bloom Filters.
 * When the current block fills up, a new block is appended with
 * capacity * CAPACITY_EXPANSION_FACTOR and target FPR * FPR_SHRINK_FACTOR.
 * Because the per-block FPRs form a geometric series, the aggregate FPR
 * stays bounded by the user-supplied P.
 *
 * Hash contract: every block delegates to bloom_filter.c, which uses
 * MurmurHash3_x64_128(key, len, seed=0) plus Normalize_Step double
 * hashing. If that contract ever changes, the on-disk format of any
 * caller that persists SBF state must be versioned; this module does
 * not persist anything itself.
 *
 * Ownership: an SBF owns its blocks and their Bloom_Filter instances.
 * SBF_Destroy frees the whole chain; SBF_Clear keeps the chain and
 * only wipes the bit arrays.
 */


#ifndef DEFAULT_EXPECTED_ELEMENT_CAPACITY

#define DEFAULT_EXPECTED_ELEMENT_CAPACITY ((uint64_t)1000000)

#endif

#ifndef CAPACITY_EXPANSION_FACTOR

#define CAPACITY_EXPANSION_FACTOR ((uint64_t)2)

#endif

#ifndef FPR_SHRINK_FACTOR

#define FPR_SHRINK_FACTOR ((double)0.5)

#endif

#ifndef MAX_BLOCK_COUNT

#define MAX_BLOCK_COUNT ((uint64_t)64)

#endif


typedef struct SBF_Block{

    Bloom_Filter* bf;

    struct SBF_Block* prev;
    struct SBF_Block* next;

}SBF_Block;

struct SBF{

    SBF_Block* head;
    SBF_Block* cur;

    uint64_t sbf_block_count;

    uint64_t sbf_ele_count;

    double sbf_tfpr;

};


//==================== Memory Management ====================

static void SBF_Block_Destroy(SBF_Block* sfbk){

    if(sfbk == NULL){
        return;
    }

    Bloom_Filter_Destroy(sfbk->bf);

    free(sfbk);

}

/*
 * SBF_Block_Create - allocate one block (descriptor + Bloom_Filter).
 *
 * On any failure, the partial state is released via SBF_Block_Destroy
 * and NULL is returned. Callers must not free the returned block
 * directly; use SBF_Block_Destroy (or SBF_Destroy for the whole chain).
 */
static SBF_Block* SBF_Block_Create(uint64_t n, double p){

    SBF_Block* sfbk = NULL;

    sfbk = (SBF_Block*)calloc(1, sizeof(SBF_Block));

    if(sfbk == NULL){
        goto cleanup;
    }

    sfbk->bf = Bloom_Filter_Create(n, p);

    if(sfbk->bf == NULL){
        goto cleanup;
    }

    return sfbk;


    cleanup:

        SBF_Block_Destroy(sfbk);
        return NULL;

}

/*
 * SBF_Create - allocate an SBF with a single head block.
 *
 * The head block's target FPR is P * (1 - FPR_SHRINK_FACTOR). For
 * FPR_SHRINK_FACTOR = 0.5 this is P/2, which is in (0, 0.5) whenever
 * P is in (0, 1); the explicit range check on initial_p is therefore
 * redundant for the default factor but stays as a guard in case a
 * caller overrides FPR_SHRINK_FACTOR via -D to a value outside [0, 1).
 *
 * On failure, both the partially-built SBF and any created block are
 * released; NULL is returned.
 */
SBF* SBF_Create(double p){

    SBF* sbf = NULL;

    if(!(p > 0.0 && p < 1.0)){
        goto cleanup;
    }

    sbf = (SBF*)calloc(1, sizeof(SBF));

    if(sbf == NULL){
        goto cleanup;
    }

    double initial_p = p * (1 - FPR_SHRINK_FACTOR);

    if(!(initial_p > 0.0 && initial_p < 1.0)){
        goto cleanup;
    }

    SBF_Block* sfbk = SBF_Block_Create(DEFAULT_EXPECTED_ELEMENT_CAPACITY, initial_p);

    if(sfbk == NULL){
        goto cleanup;
    }

    sbf->head            = sfbk;
    sbf->cur             = sfbk;
    sbf->sbf_block_count = 1;
    sbf->sbf_tfpr        = p;

    return sbf;


    cleanup:

        SBF_Destroy(sbf);
        return NULL;

}

/*
 * SBF_Expand - append a new block with larger capacity and lower FPR.
 *
 * Returns false (and leaves the chain unchanged) in any of these cases:
 *   1. the chain already has MAX_BLOCK_COUNT blocks;
 *   2. capacity * CAPACITY_EXPANSION_FACTOR would overflow uint64_t;
 *   3. the new target FPR would leave the open interval (0, 1);
 *   4. allocating the new block fails (including Bloom_Filter_Create
 *      failing because the recomputed bit size is too large).
 *
 * On success, the new block becomes the tail and sbf->cur; the old
 * tail's next pointer is updated. The caller (SBF_Element_Mark) is
 * responsible for retrying the insert after a successful expand.
 */
static bool SBF_Expand(SBF* sbf){

    if(sbf == NULL || sbf->head == NULL || sbf->cur == NULL){
        return false;
    }

    if(sbf->sbf_block_count >= MAX_BLOCK_COUNT){
        return false;
    }

    uint64_t prev_capacity = Bloom_Filter_Capacity(sbf->cur->bf);

    if(prev_capacity > UINT64_MAX / CAPACITY_EXPANSION_FACTOR){
        return false;
    }

    uint64_t new_capacity = prev_capacity * CAPACITY_EXPANSION_FACTOR;
    double new_p          = Bloom_Filter_TFPR(sbf->cur->bf) * FPR_SHRINK_FACTOR;

    if(!(new_p > 0.0 && new_p < 1.0)){
        return false;
    }

    SBF_Block* sfbk = SBF_Block_Create(new_capacity, new_p);

    if(sfbk == NULL){
        return false;
    }

    sfbk->prev = sbf->cur;
    sbf->cur->next = sfbk;
    sbf->cur = sfbk;

    sbf->sbf_block_count++;

    return true;

}

/*
 * SBF_Clear - wipe every block's bit array and reset the global
 * element counter.
 *
 * This does NOT free blocks, does NOT change sbf_block_count, and
 * does NOT change sbf_tfpr: they are design constants / structural
 * state. sbf->cur is reset to sbf->head, so the next insert goes
 * back into the first block. Use SBF_Destroy to fully release the
 * chain.
 */
void SBF_Clear(SBF* sbf){

    if(sbf == NULL){
        return;
    }

    SBF_Block* sfbk = sbf->head;

    while(sfbk != NULL){

        Bloom_Filter_Clear(sfbk->bf);
        sfbk = sfbk->next;

    }

    sbf->cur = sbf->head;
    sbf->sbf_ele_count = 0;

}

/*
 * SBF_Destroy - free every block and the SBF handle itself.
 *
 * Walks the chain via next from head and releases each block with
 * SBF_Block_Destroy. prev pointers are not used here; they exist
 * only to support SBF_Element_Query. Passing NULL is a silent no-op.
 */
void SBF_Destroy(SBF* sbf){

    if(sbf == NULL){
        return;
    }

    while(sbf->head != NULL){

        SBF_Block* sfbk = sbf->head;
        sbf->head = sbf->head->next;

        SBF_Block_Destroy(sfbk);

    }

    free(sbf);

}


//==================== Operating Functions ====================

/*
 * SBF_Element_Mark - insert one byte sequence into the current block.
 *
 * Returns false if:
 *   - sbf is NULL, or
 *   - key is NULL and len > 0 (invalid; key == NULL with len == 0 is
 *     a valid "empty key", same as bloom_filter.c), or
 *   - the current block is full and SBF_Expand fails (see the four
 *     failure modes documented on SBF_Expand).
 *
 * Note: a false return does NOT mean "element already present".
 * Membership must be tested with SBF_Element_Query. The element
 * counter is incremented only when the underlying Bloom_Filter_Mark
 * succeeds, so a failed expansion does not silently bump the count.
 */
bool SBF_Element_Mark(SBF* sbf, const void* key, uint64_t len){

    if(sbf == NULL || (key == NULL && len > 0)){
        return false;
    }

    if(Bloom_Filter_IsFull(sbf->cur->bf) && SBF_Expand(sbf) == false){
        return false;
    }

    if(Bloom_Filter_Element_Mark(sbf->cur->bf, key, len)){
        sbf->sbf_ele_count++;
        return true;
    }

    return false;

}

/*
 * SBF_Element_Query - test whether a byte sequence is in any block.
 *
 * Walks from sbf->cur (newest block) backwards through prev to head
 * (oldest block). Newest blocks have lower FPR, so a hit is most
 * likely there; returning on the first hit avoids touching older,
 * more saturated blocks. The query is read-only with respect to the
 * bit arrays but is not safe concurrently with a writer (see the
 * header's thread-safety note).
 *
 * Same key/len semantics as SBF_Element_Mark: key == NULL with
 * len == 0 is a valid empty key.
 */
bool SBF_Element_Query(const SBF* const sbf, const void* key, uint64_t len){

    if(sbf == NULL || (key == NULL && len > 0)){
        return false;
    }

    const SBF_Block* sfbk = sbf->cur;

    while(sfbk != NULL){

        if(Bloom_Filter_Element_Query(sfbk->bf, key, len)){
            return true;
        }

        sfbk = sfbk->prev;

    }

    return false;

}


//==================== Parameter Accessor Functions ====================

uint64_t SBF_Block_Count(const SBF* const sbf){

    if(sbf == NULL){
        return UINT64_MAX;
    }

    return sbf->sbf_block_count;

}

uint64_t SBF_Ele_Count(const SBF* const sbf){

    if(sbf == NULL){
        return UINT64_MAX;
    }

    return sbf->sbf_ele_count;

}

/*
 * SBF_CFPR - estimate the current aggregate false-positive rate.
 *
 * Each block contributes its *current* p_i, computed by
 * Bloom_Filter_CFPR(block) from that block's live ele_count,
 * hash_count, and bit_size. The aggregate is
 *
 *     P_current = 1 - prod_i (1 - p_i)
 *
 * which is the probability that a non-member query hits at least one
 * block. This is an *estimate*, not a measurement: it uses the
 * standard Bloom-Filter FPR formula, which assumes the hash functions
 * are independent and uniform. The measured rate on a finite sample
 * may deviate due to sampling noise, but it is not an "upper bound"
 * in any formal sense.
 *
 * Numerical note: the product is evaluated in log space
 * (-expm1(sum log1p(-p_i))) to avoid underflow and cancellation when
 * many blocks accumulate.
 *
 * p_i == 1.0 is handled explicitly: log1p(-1) would be -inf, and the
 * whole chain would be saturated, so CFPR returns 1.0.
 */
double SBF_CFPR(const SBF* const sbf){

    if(sbf == NULL){
        return NAN;
    }

    double log_miss_sum = 0.0;

    const SBF_Block* sfbk = sbf->head;

    while(sfbk != NULL){

        double p_i = Bloom_Filter_CFPR(sfbk->bf);

        if(!(p_i >= 0.0 && p_i <= 1.0)){
            return NAN;
        }

        /* p_i == 1.0 would make log1p(-1) = -inf; the whole filter would
         * then be saturated (CFPR = 1). Handle it explicitly. */
        if(p_i == 1.0){
            return 1.0;
        }

        log_miss_sum += log1p(-p_i);

        sfbk = sfbk->next;

    }

    /* log_miss_sum <= 0 always (each term is <= 0). */
    double p_total = -expm1(log_miss_sum);

    /* Clamp defensively: rounding may push the value marginally outside [0,1]. */
    if(p_total < 0.0){
        p_total = 0.0; 
    }

    if(p_total > 1.0){
        p_total = 1.0;
    }

    return p_total;

}

double SBF_TFPR(const SBF* const sbf){

    if(sbf == NULL){
        return NAN;
    }

    return sbf->sbf_tfpr;

}