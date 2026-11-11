/**
 * @file bloom_filter.h
 * @brief Public interface of a classic Bloom Filter.
 *
 * A Bloom Filter is a space-efficient probabilistic data structure used to
 * test **whether an element may be present in a set**. It has two key
 * properties:
 *
 * - **No false negatives**: querying an inserted element always returns `true`.
 * - **Possible false positives**: querying a non-inserted element may return
 *   `true`, with probability bounded by the target false-positive rate @c p
 *   specified at creation time.
 *
 * This implementation uses a bit array plus double hashing to derive k hash
 * positions: the i-th position is computed as `h1 + i * h2 (mod m)`, avoiding
 * the need to store k independent hash functions.
 *
 * @par Typical usage
 * @code
 * Bloom_Filter* bf = Bloom_Filter_Create(1000000, 0.01);
 * if (bf == NULL) { ... }
 *
 * Bloom_Filter_Element_Mark(bf, "user_42", 7);
 * bool maybe = Bloom_Filter_Element_Query(bf, "user_42", 7);  // always true
 *
 * Bloom_Filter_Destroy(bf);
 * @endcode
 *
 * @par Thread safety
 * This module is **not thread-safe**. Callers must provide external locking
 * in multithreaded environments.
 *
 * @warning A Bloom Filter does not support deleting individual elements
 *          (`Clear` wipes the entire filter).
 *
 * @author  (to be filled)
 * @version 1.0
 * @date    2026.10.28
 */
#pragma once

#include <stdbool.h>    /**< bool type */
#include <stdint.h>     /**< uint64_t type */

/**
 * @brief Opaque handle to a Bloom Filter.
 *
 * The actual fields are defined in `bloom_filter.c`; callers neither need nor
 * are able to access the internal structure.
 */
typedef struct Bloom_Filter Bloom_Filter;


/* ==================================================================
 *  Memory Management
 * ================================================================== */

/**
 * @brief Creates a Bloom Filter.
 *
 * Given the expected number of elements @p n and the target false-positive
 * rate @p p, the optimal bit-array size @c m and hash count @c k are derived
 * automatically:
 * @f[
 *     m = -\frac{n \ln p}{(\ln 2)^2}, \qquad
 *     k = -\frac{\ln p}{\ln 2}
 * @f]
 *
 * @param[in] n Expected number of elements to insert; must be @c > 0.
 * @param[in] p Target false-positive rate; must satisfy @c 0.0 < p < 1.0.
 *
 * @return A newly created handle on success; @c NULL on invalid arguments or
 *         allocation failure.
 *
 * @retval NULL @p n == 0
 * @retval NULL @p n > UINT64_MAX / 2 (the element counter could wrap)
 * @retval NULL @p p is outside the open interval (0, 1)
 * @retval NULL Out of memory, or the bit-array size exceeds the platform's
 *              representable range
 *
 * @note The returned handle must be released with Bloom_Filter_Destroy(),
 *       otherwise the memory leaks.
 * @see Bloom_Filter_Destroy, Bloom_Filter_Clear
 */
Bloom_Filter* Bloom_Filter_Create(uint64_t n, double p);

/**
 * @brief Clears all data stored in the filter.
 *
 * Sets every bit in the bit array to 0 and resets the element counter to 0.
 * Metadata such as capacity and target false-positive rate are **preserved**;
 * the filter can be reused immediately after clearing.
 *
 * @param[in,out] bf Filter handle. If @c NULL, this is a silent no-op.
 *
 * @post Bloom_Filter_Element_Count(bf) == 0
 * @post Bloom_Filter_IsFull(bf) == false
 *
 * @see Bloom_Filter_Destroy
 */
void Bloom_Filter_Clear(Bloom_Filter* bf);

/**
 * @brief Destroys the filter and releases all associated memory.
 *
 * @param[in] bf Handle to destroy. If @c NULL, this is a silent no-op.
 *
 * @warning After this call @p bf is invalid; any subsequent access is
 *          undefined behavior (dangling pointer).
 * @see Bloom_Filter_Create
 */
void Bloom_Filter_Destroy(Bloom_Filter* bf);


/* ==================================================================
 *  Operations
 * ================================================================== */

/**
 * @brief Inserts a byte sequence into the filter.
 *
 * Computes k hash positions for @p key, sets the corresponding bits to 1, and
 * increments the element counter by one.
 *
 * @param[in,out] bf  Filter handle.
 * @param[in]     key Pointer to the bytes. May be @c NULL only when @p len == 0.
 * @param[in]     len Number of bytes to insert.
 *
 * @note `key` is treated as a raw byte sequence of `len` bytes; it does
 *       NOT need to be NUL-terminated. Any binary data is allowed, including
 *       embedded '\0' bytes.
 *
 * @return @c true on success; @c false on invalid arguments.
 *
 * @retval false @p bf == NULL
 * @retval false @p key == NULL and @p len > 0
 *
 * @note Re-inserting the same byte sequence is **idempotent at the bit level**, but
 *       it **does increment the element counter again** (the counter semantics
 *       are "number of insert calls", not "number of distinct elements").
 * @note Even when the filter has reached its capacity
 *       (Bloom_Filter_IsFull returns @c true), this function still succeeds
 *       and increments the counter. This is the documented behavior of a
 *       classic Bloom Filter: capacity is a statistical design target, not a
 *       hard upper bound.
 *
 * @see Bloom_Filter_Element_Query, Bloom_Filter_IsFull
 */
bool Bloom_Filter_Element_Mark(Bloom_Filter* bf, const void* key, uint64_t len);

/**
 * @brief Queries whether a byte sequence **may be present** in the filter.
 *
 * @param[in] bf  Filter handle.
 * @param[in] key Pointer to the bytes. May be @c NULL only when @p len == 0.
 * @param[in] len Number of bytes to insert.
 *
 * @note `key` is treated as a raw byte sequence of `len` bytes; it does
 *       NOT need to be NUL-terminated. Any binary data is allowed, including
 *       embedded '\0' bytes.
 *
 * @return The query result.
 * @retval true  The byte sequence **may** be present (false positives are possible).
 * @retval false The byte sequence is **definitely not** present (no false negatives).
 * @retval false Invalid arguments (@p bf or @p key == NULL and @p len > 0).
 *
 * @note A @c false return means "definitely absent" and can be safely used
 *       for short-circuit filtering.
 * @note A @c true return does not imply actual presence; the caller must
 *       verify with an authoritative source when needed.
 *
 * @see Bloom_Filter_Element_Mark
 */
bool Bloom_Filter_Element_Query(const Bloom_Filter* const bf, const void* key, uint64_t len);


/* ==================================================================
 *  Parameter Accessors
 * ================================================================== */

/**
 * @brief Returns the cumulative count of insertions.
 *
 * @param[in] bf Filter handle.
 * @return Cumulative insertion count; @c UINT64_MAX as an error sentinel if
 *         @p bf == NULL.
 *
 * @warning The count includes duplicate insertions and may exceed the value
 *          returned by Bloom_Filter_Capacity().
 * @see Bloom_Filter_Capacity
 */
uint64_t Bloom_Filter_Element_Count(const Bloom_Filter* const bf);

/**
 * @brief Returns the expected element capacity @p n specified at creation.
 *
 * @param[in] bf Filter handle.
 * @return Expected capacity; @c UINT64_MAX as an error sentinel if
 *         @p bf == NULL.
 */
uint64_t Bloom_Filter_Capacity(const Bloom_Filter* const bf);

/**
 * @brief Returns the **target** false-positive rate (Target FPR).
 *
 * This equals the @p p argument passed to Bloom_Filter_Create(). It is a
 * design constant and does not change as elements are inserted.
 *
 * @param[in] bf Filter handle.
 * @return Target false-positive rate; @c NAN if @p bf == NULL.
 *
 * @see Bloom_Filter_CFPR
 */
double Bloom_Filter_TFPR(const Bloom_Filter* const bf);

/**
 * @brief Computes the theoretical estimate of the **current** false-positive
 *        rate (Current FPR).
 *
 * Uses the formula:
 * @f[
 *     p_{\text{current}} = \left(1 - e^{-kn/m}\right)^{k}
 * @f]
 * where @c k is the number of hash functions, @c n is the current element
 * count, and @c m is the bit-array size.
 *
 * @param[in] bf Filter handle.
 * @return Estimated current false-positive rate in [0, 1]; @c NAN if
 *         @p bf == NULL.
 *
 * @note This value is a **theoretical estimate**; measured rates may deviate
 *       due to sampling noise.
 * @note An empty filter returns 0.
 * @see Bloom_Filter_TFPR
 */
double Bloom_Filter_CFPR(const Bloom_Filter* const bf);

/**
 * @brief Checks whether the filter has reached its design capacity.
 *
 * Returns @c true when Bloom_Filter_Element_Count() >= Bloom_Filter_Capacity().
 *
 * @param[in] bf Filter handle.
 * @return Whether the filter is full.
 * @retval true  Full, or @p bf == NULL (conservative: prevents misuse of a
 *               null handle).
 * @retval false Capacity not yet reached.
 *
 * @note After reaching capacity, Bloom_Filter_Element_Mark() **still works**,
 *       but the actual false-positive rate will exceed the design target.
 * @see Bloom_Filter_Element_Mark
 */
bool Bloom_Filter_IsFull(const Bloom_Filter* const bf);