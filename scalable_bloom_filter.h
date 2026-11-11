/**
 * @file scalable_bloom_filter.h
 * @brief Public interface of a Scalable Bloom Filter (SBF).
 *
 * A classic Bloom Filter requires its capacity to be fixed at creation time;
 * once the number of elements approaches that capacity, the false-positive
 * rate rises sharply. A **Scalable Bloom Filter** overcomes this limitation
 * by maintaining a **block chain**: when the current block fills up, a new
 * block is appended automatically:
 *
 * - New-block capacity grows by a constant factor
 *   @c CAPACITY_EXPANSION_FACTOR (default: 2×).
 * - New-block target false-positive rate shrinks by
 *   @c FPR_SHRINK_FACTOR (default: ×0.5), keeping the **aggregate
 *   false-positive rate** of the whole chain bounded by the user-specified
 *   @p p.
 * - Queries traverse the chain from tail (newest block) to head (oldest
 *   block), returning @c true on the first hit.
 *
 * @par Aggregate false-positive rate
 * Let @f$ p_i @f$ be the target false-positive rate of block i. The aggregate
 * false-positive rate of the chain is
 * @f[
 *     P_{\text{total}} = 1 - \prod_{i=0}^{N-1} \left(1 - p_i\right)
 * @f]
 * Because @f$ p_i @f$ forms a geometric sequence, the series converges and the
 * aggregate error rate stays bounded as the chain grows.
 *
 * @par Typical usage
 * @code
 * SBF* sbf = SBF_Create(0.01);          // aggregate target FPR = 1%
 * if (sbf == NULL) { ... }
 *
 * SBF_Element_Mark(sbf, "user_42", 7);     // auto-expands as needed
 * bool maybe = SBF_Element_Query(sbf, "user_42", 7);
 *
 * SBF_Destroy(sbf);
 * @endcode
 *
 * @par Thread safety
 * This module is **not thread-safe**. Expansion mutates the chain structure,
 * so concurrent readers/writers must be externally locked.
 * 
 * @par Hash contract
 * Every block delegates to bloom_filter.c, which uses
 * MurmurHash3_x64_128(key, len, seed=0) plus Normalize_Step double
 * hashing. If a caller persists SBF state across process restarts,
 * the persisted format must be versioned together with this contract;
 * this module does not persist anything itself.
 *
 * @par Key semantics
 * `key` / `len` follow bloom_filter.h: `key` is a raw byte sequence
 * of `len` bytes, not NUL-terminated; embedded '\0' bytes are allowed.
 * `key == NULL` is valid only when `len == 0` (empty key).
 *
 * @par Expansion is not rollback-safe
 * If SBF_Element_Mark triggers an expansion and the subsequent insert
 * fails, the newly appended (empty) block stays in the chain. This is
 * harmless for correctness (queries will simply find nothing there)
 * but means block_count can grow even when Mark returns false.
 *
 * @note This module depends on the classic Bloom Filter in `bloom_filter.h`
 *       for its per-block implementation.
 *
 * @author  (to be filled)
 * @version 1.0
 * @date    2026.10.28
 */
#pragma once

#include "bloom_filter.h"

#include <stdbool.h>    /**< bool type */


/**
 * @brief Opaque handle to a Scalable Bloom Filter.
 *
 * Internally structured as a doubly linked list of SBF_Block nodes. Fields are
 * defined in `scalable_bloom_filter.c`; callers do not need access to the
 * internal structure.
 */
typedef struct SBF SBF;


/* ==================================================================
 *  Memory Management
 * ================================================================== */

/**
 * @brief Creates a Scalable Bloom Filter.
 *
 * Initializes a "head block" with capacity @c DEFAULT_EXPECTED_ELEMENT_CAPACITY
 * (default: 1,000,000) and a target false-positive rate of
 * @f$ p \times (1 - \text{FPR\_SHRINK\_FACTOR}) @f$.
 *
 * @param[in] p **Aggregate** target FPR for the whole chain. This is
 *              not the per-block FPR: the head block gets p/2, the
 *              next p/4, and so on, so the union bound equals p.
 *              Must satisfy @c 0.0 < p < 1.0.
 *
 * @return A newly created SBF handle on success; @c NULL on invalid arguments
 *         or allocation failure.
 *
 * @retval NULL @p p is outside the open interval (0, 1)
 * @retval NULL Out of memory
 *
 * @note The returned handle must be released with SBF_Destroy().
 * @see SBF_Destroy, SBF_Clear
 */
SBF* SBF_Create(double p);

/**
 * @brief Clears all data in the chain but **keeps the allocated block
 *        structure**.
 *
 * Iterates over the entire chain, calling Bloom_Filter_Clear() on each block;
 * then resets the write cursor to the head block and zeroes the global element
 * counter.
 *
 * @param[in,out] sbf SBF handle. If @c NULL, this is a silent no-op.
 *
 * @post SBF_Ele_Count(sbf) == 0
 * @post SBF_Block_Count(sbf) is **unchanged** (blocks are not freed)
 * @post The next insertion writes into the head block
 * @post SBF_TFPR(sbf) is unchanged (it is a design constant)
 *
 * @note Blocks are retained to avoid performance jitter from repeated
 *       re-expansion. Use SBF_Destroy() to fully release memory.
 * @see SBF_Destroy
 */
void SBF_Clear(SBF* sbf);

/**
 * @brief Destroys the filter, freeing all blocks and their bit arrays.
 *
 * @param[in] sbf Handle to destroy. If @c NULL, this is a silent no-op.
 *
 * @warning After this call @p sbf is invalid; any subsequent access is
 *          undefined behavior (dangling pointer).
 * @post SBF_TFPR(sbf) is unchanged (it is a design constant)
 * 
 * @see SBF_Create
 */
void SBF_Destroy(SBF* sbf);


/* ==================================================================
 *  Operations
 * ================================================================== */

/**
 * @brief Inserts a byte sequence into the scalable filter.
 *
 * If the current write block (@c cur) is full, an automatic expansion is
 * triggered first (appending a new block with larger capacity and a lower
 * false-positive rate), after which the element is written into the new
 * block.
 *
 * @param[in,out] sbf SBF handle.
 * @param[in]     key Pointer to the bytes. May be @c NULL only when
 *                    @p len == 0.
 * @param[in]     len Number of bytes to insert.
 *
 * @return @c true on success; @c false on invalid arguments or expansion
 *         failure.
 *
 * @retval false @p sbf == NULL
 * @retval false @p key == NULL and @p len > 0
 * @retval false Expansion failed (MAX_BLOCK_COUNT reached, capacity
 *               overflow, new FPR outside (0, 1), or allocation failed)
 *
 * @note Insertions always target the **current write block** (the chain tail);
 *       historical blocks are no longer written to.
 * @see SBF_Element_Query
 */
bool SBF_Element_Mark(SBF* sbf, const void* key, uint64_t len);

/**
 * @brief Queries whether a byte sequence **may be present** anywhere in the chain.
 *
 * Traverses the chain from tail (newest block) to head (oldest block),
 * calling Bloom_Filter_Element_Query() on each block. Returns @c true on the
 * first hit. Because newer blocks have lower false-positive rates, starting
 * from the tail allows **recently inserted elements to be found fastest**.
 *
 * @param[in] sbf SBF handle (const-qualified; internal state is not
 *                modified).
 * @param[in] key Pointer to the bytes. May be @c NULL only when
 *                @p len == 0.
 * @param[in] len Number of bytes to query.
 *
 * @return The query result.
 * @retval true  The byte sequence **may** exist in some block of the chain (false
 *               positives are possible).
 * @retval false @p sbf == NULL
 * @retval false @p key == NULL and @p len > 0
 *
 * @see SBF_Element_Mark
 */
bool SBF_Element_Query(const SBF* const sbf, const void* key, uint64_t len);


/* ==================================================================
 *  Parameter Accessors
 * ================================================================== */

/**
 * @brief Returns the current number of blocks in the chain.
 *
 * @param[in] sbf SBF handle.
 * @return Block count; @c UINT64_MAX as an error sentinel if @p sbf == NULL.
 *
 * @note The block count grows monotonically with expansions and only resets
 *       to zero after SBF_Destroy() (at which point the handle is invalid).
 * @note Block count grows on every successful SBF_Expand. It can
 *       also grow when SBF_Element_Mark triggers an expansion and
 *       the subsequent insert fails (the empty block stays in the
 *       chain); see the "Expansion is not rollback-safe" note above.
 */
uint64_t SBF_Block_Count(const SBF* const sbf);

/**
 * @brief Returns the cumulative number of elements inserted across the chain.
 *
 * @param[in] sbf SBF handle.
 * @return Cumulative insertion count; @c UINT64_MAX as an error sentinel if
 *         @p sbf == NULL.
 *
 * @note The count includes duplicate insertions and may exceed the initial
 *       block capacity.
 * @note Block count grows on every successful SBF_Expand. It can also
 *       grow when SBF_Element_Mark triggers an expansion and the
 *       subsequent insert fails (the empty block stays in the chain).
 * @see SBF_Block_Count
 */
uint64_t SBF_Ele_Count(const SBF* const sbf);

/**
 * @brief Computes the **current** aggregate false-positive rate estimate of
 *        the whole chain.
 *
 * Uses the formula:
 * @f[
 *     P_{\text{current}} = 1 - \prod_{i=0}^{N-1} \left(1 - p_i^{\text{cur}}\right)
 * @f]
 * where @f$ p_i^{\text{cur}} @f$ is the current false-positive rate estimate
 * of block i as returned by Bloom_Filter_CFPR().
 *
 * @param[in] sbf SBF handle.
 * @return Aggregate current false-positive rate in [0, 1]; @c NAN if
 *         @p sbf == NULL or any block returns an invalid value.
 *
 * @note An empty filter (all blocks empty) returns 0.
 * @see SBF_TFPR
 */
double SBF_CFPR(const SBF* const sbf);

/**
 * @brief Returns the **aggregate target** false-positive rate (Target FPR)
 *        specified at creation.
 *
 * This equals the @p p argument passed to SBF_Create(). It is a design
 * constant for the whole chain and does not change as elements are inserted or
 * the chain expands.
 *
 * @param[in] sbf SBF handle.
 * @return Aggregate target false-positive rate; @c NAN if @p sbf == NULL.
 *
 * @see SBF_CFPR
 */
double SBF_TFPR(const SBF* const sbf);