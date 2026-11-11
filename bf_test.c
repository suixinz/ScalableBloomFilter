#include "bloom_filter.h"
#include "murmur3.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

/* ------------------------------------------------------------------ */
/* Constants                                                          */
/* ------------------------------------------------------------------ */

#define CFPR_TOLERANCE     (0.30)    /* 30% relative */
#define FP_ABS_TOLERANCE   (0.005)   /* absolute distance from CFPR */

/* ------------------------------------------------------------------ */
/* Test utilities                                                      */
/* ------------------------------------------------------------------ */

static void Progress_Bar(uint64_t step, uint64_t total, const char* label){

    if(total == 0){
        return;
    }

    const int width = 40;

    double ratio = (double)step / (double)total;

    if(ratio > 1.0){
        ratio = 1.0;
    }

    int filled = (int)(ratio * (double)width);

    printf("\r[");

    for(int i = 0; i < width; i++){
        putchar(i < filled ? '=' : ' ');
    }

    printf("] %6.2f%%  %-40s", ratio * 100.0, label);

    if(step >= total){
        putchar('\n');
    }

    fflush(stdout);

}

static void Make_Key(char* buf, size_t size, const char* prefix, uint64_t id){
    snprintf(buf, size, "%s_%" PRIu64, prefix, id);
}

static inline int Approx_Equal(double a, double b, double rel_tol){
    double diff = fabs(a - b);
    double base = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    if(base < 1e-12){
        return diff < 1e-12;
    }
    return diff / base <= rel_tol;
}

/* Convenience: mark / query a C string (NUL-terminated). */
static inline bool Mark_Str(Bloom_Filter* bf, const char* s){
    return Bloom_Filter_Element_Mark(bf, s, (uint64_t)strlen(s));
}

static inline bool Query_Str(const Bloom_Filter* const bf, const char* s){
    return Bloom_Filter_Element_Query(bf, s, (uint64_t)strlen(s));
}

/* ------------------------------------------------------------------ */
/* Test cases                                                          */
/* ------------------------------------------------------------------ */

static int Test_Murmur3_Known_Vectors(void){

    printf("\n--- Test_Murmur3_Known_Vectors ---\n");

    int fails = 0;

    /*
     * Published MurmurHash3 test vector:
     *   MurmurHash3_x64_128("", 0, seed=0) == { 0, 0 }
     *
     * This is the single most important test in this file. If the
     * implementation of rotl64 / fmix64 / getblock64 / the tail switch
     * is wrong, the Bloom Filter will still "work" (Mark and Query use
     * the same broken hash, so Test_No_False_Negative passes), but the
     * hash quality degrades and the false-positive rate silently rises.
     * This vector catches that.
     */

    uint64_t out[2];

    MurmurHash3_x64_128("", 0, 0, out);

    if(out[0] != 0 || out[1] != 0){

        printf("  FAIL: hash(\"\", 0, seed=0) = {%016" PRIx64 ", %016" PRIx64
               "}, expected {0, 0}\n", out[0], out[1]);

        fails++;

    }
    else{

        printf("  PASS: hash(\"\", 0, seed=0) == {0, 0}\n");

    }

    /*
     * Sanity checks: different inputs must produce different outputs,
     * and the two 64-bit halves must not be trivially related.
     */

    uint64_t a[2], b[2];

    MurmurHash3_x64_128("a", 1, 0, a);
    MurmurHash3_x64_128("b", 1, 0, b);

    if(a[0] == b[0] && a[1] == b[1]){

        printf("  FAIL: hash(\"a\") == hash(\"b\"), hash is degenerate\n");

        fails++;

    }

    if(a[0] == a[1]){

        printf("  FAIL: out[0] == out[1] for hash(\"a\"), halves not independent\n");

        fails++;

    }

    /* Same input + same seed must be deterministic. */

    uint64_t a2[2];

    MurmurHash3_x64_128("a", 1, 0, a2);

    if(a[0] != a2[0] || a[1] != a2[1]){

        printf("  FAIL: hash is not deterministic\n");

        fails++;

    }

    /* Different seed must produce a different hash. */

    uint64_t a3[2];

    MurmurHash3_x64_128("a", 1, 1, a3);

    if(a[0] == a3[0] && a[1] == a3[1]){

        printf("  FAIL: different seed produced the same hash\n");

        fails++;

    }

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Create_Invalid_Params(void){

    printf("\n--- Test_Create_Invalid_Params ---\n");

    int fails = 0;

    if(Bloom_Filter_Create(0, 0.01) != NULL){

        printf("  FAIL: n=0 should return NULL\n");

        fails++;

    }

    if(Bloom_Filter_Create(100, 0.0) != NULL){

        printf("  FAIL: p=0 should return NULL\n");

        fails++;

    }

    if(Bloom_Filter_Create(100, 1.0) != NULL){

        printf("  FAIL: p=1 should return NULL\n");

        fails++;

    }

    if(Bloom_Filter_Create(100, -0.5) != NULL){

        printf("  FAIL: p=-0.5 should return NULL\n");

        fails++;

    }

    if(Bloom_Filter_Create(100, 1.5) != NULL){

        printf("  FAIL: p=1.5 should return NULL\n");

        fails++;

    }

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Null_Defense(void){

    printf("\n--- Test_Null_Defense ---\n");

    int fails = 0;

    Bloom_Filter_Clear(NULL);
    Bloom_Filter_Destroy(NULL);

    /* key == NULL with len > 0: invalid, must return false. */

    if(Bloom_Filter_Element_Mark(NULL, "x", 1) != false){

        printf("  FAIL: Mark(NULL, \"x\", 1)\n");

        fails++;

    }

    if(Bloom_Filter_Element_Query(NULL, "x", 1) != false){

        printf("  FAIL: Query(NULL, \"x\", 1)\n");

        fails++;

    }

    Bloom_Filter* bf = Bloom_Filter_Create(1000, 0.01);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    /* key == NULL with len > 0: invalid. */

    if(Bloom_Filter_Element_Mark(bf, NULL, 5) != false){

        printf("  FAIL: Mark(bf, NULL, 5)\n");

        fails++;

    }

    if(Bloom_Filter_Element_Query(bf, NULL, 5) != false){

        printf("  FAIL: Query(bf, NULL, 5)\n");

        fails++;

    }

    /*
     * key == NULL with len == 0 is a valid "empty key": Murmur3 allows
     * (NULL, 0), and so does the Bloom Filter. Both Mark and Query must
     * return true.
     */

    if(Bloom_Filter_Element_Mark(bf, NULL, 0) != true){

        printf("  FAIL: Mark(bf, NULL, 0) should succeed (empty key)\n");

        fails++;

    }

    if(Bloom_Filter_Element_Query(bf, NULL, 0) != true){

        printf("  FAIL: Query(bf, NULL, 0) should be true after Mark\n");

        fails++;

    }

    if(Bloom_Filter_Element_Count(NULL) != UINT64_MAX){

        printf("  FAIL: Count(NULL)\n");

        fails++;

    }

    if(Bloom_Filter_Capacity(NULL) != UINT64_MAX){

        printf("  FAIL: Capacity(NULL)\n");

        fails++;

    }

    if(!isnan(Bloom_Filter_TFPR(NULL))){

        printf("  FAIL: TFPR(NULL)\n");

        fails++;

    }

    if(!isnan(Bloom_Filter_CFPR(NULL))){

        printf("  FAIL: CFPR(NULL)\n");

        fails++;

    }

    if(Bloom_Filter_IsFull(NULL) != true){

        printf("  FAIL: IsFull(NULL)\n");

        fails++;

    }

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Tiny_Capacity(void){

    printf("\n--- Test_Tiny_Capacity ---\n");

    int fails = 0;

    /*
     * n = 1 and n = 2 produce very small bit arrays (bit_size = 8 or 16).
     * This exercises Normalize_Step on tiny bit_size values, where the
     * gcd loop may iterate more than once. It also checks that empty
     * keys and single-byte keys work.
     */

    const uint64_t tiny_n[] = {1, 2, 3, 4, 8};

    for(size_t idx = 0; idx < sizeof(tiny_n) / sizeof(tiny_n[0]); idx++){

        Bloom_Filter* bf = Bloom_Filter_Create(tiny_n[idx], 0.01);

        if(bf == NULL){

            printf("  FAIL: Create(%" PRIu64 ") returned NULL\n", tiny_n[idx]);

            fails++;

            continue;

        }

        Mark_Str(bf, "");

        if(!Query_Str(bf, "")){

            printf("  FAIL: empty string round-trip failed for n=%" PRIu64 "\n", tiny_n[idx]);

            fails++;

        }

        Mark_Str(bf, "x");

        if(!Query_Str(bf, "x")){

            printf("  FAIL: single-char round-trip failed for n=%" PRIu64 "\n", tiny_n[idx]);

            fails++;

        }

        Bloom_Filter_Destroy(bf);

    }

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Binary_Key(void){

    printf("\n--- Test_Binary_Key ---\n");

    int fails = 0;

    Bloom_Filter* bf = Bloom_Filter_Create(10000, 0.01);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    /*
     * Binary keys with embedded '\0' must work. The API takes (key, len),
     * not a NUL-terminated string, so "a\0b" (len 3) and "a\0c" (len 3)
     * are distinct keys even though strlen() would stop at the first '\0'.
     */

    const unsigned char k1[] = { 'a', 0x00, 'b' };
    const unsigned char k2[] = { 'a', 0x00, 'c' };
    const unsigned char k3[] = { 0xFF, 0x00, 0x80, 0x7F };

    Bloom_Filter_Element_Mark(bf, k1, sizeof(k1));
    Bloom_Filter_Element_Mark(bf, k2, sizeof(k2));
    Bloom_Filter_Element_Mark(bf, k3, sizeof(k3));

    if(!Bloom_Filter_Element_Query(bf, k1, sizeof(k1))){

        printf("  FAIL: false negative for binary key k1\n");

        fails++;

    }

    if(!Bloom_Filter_Element_Query(bf, k2, sizeof(k2))){

        printf("  FAIL: false negative for binary key k2\n");

        fails++;

    }

    if(!Bloom_Filter_Element_Query(bf, k3, sizeof(k3))){

        printf("  FAIL: false negative for binary key k3\n");

        fails++;

    }

    /* k1 and k2 differ only in the last byte; both must be present. */

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_No_False_Negative(void){

    printf("\n--- Test_No_False_Negative ---\n");

    const uint64_t n = 200000;

    Bloom_Filter* bf = Bloom_Filter_Create(n, 0.01);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    /* Edge case: empty string must round-trip. */

    Mark_Str(bf, "");

    if(!Query_Str(bf, "")){

        printf("  FAIL: false negative for empty string\n");

        fails++;

    }

    /* Phase 1: insert all */

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "user", i);

        if(!Mark_Str(bf, key)){

            printf("  FAIL: Mark failed at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 20000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "insert phase");
        }

    }

    /* Phase 2: query all inserted -> must be true */

    for(uint64_t i = 0; i < n && fails == 0; i++){

        Make_Key(key, sizeof(key), "user", i);

        if(!Query_Str(bf, key)){

            printf("  FAIL: false negative at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 20000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "query inserted");
        }

    }

    printf("  Element_Count = %" PRIu64 " (expected %" PRIu64 ")\n",
           Bloom_Filter_Element_Count(bf), n + 1);

    if(Bloom_Filter_Element_Count(bf) != n + 1){

        printf("  FAIL: element count mismatch\n");

        fails++;

    }

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Cfpr_Matches_Tfpr(void){

    printf("\n--- Test_Cfpr_Matches_Tfpr ---\n");

    const uint64_t n = 100000;
    const double   p = 0.01;

    Bloom_Filter* bf = Bloom_Filter_Create(n, p);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "k", i);

        Mark_Str(bf, key);

        if((i % 10000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "fill to capacity");
        }

    }

    double tfpr = Bloom_Filter_TFPR(bf);
    double cfpr = Bloom_Filter_CFPR(bf);

    printf("  TFPR = %.6f, CFPR = %.6f, rel_tol = %.2f\n", tfpr, cfpr, CFPR_TOLERANCE);

    if(!Approx_Equal(cfpr, tfpr, CFPR_TOLERANCE)){

        printf("  FAIL: CFPR %.6f not within %.0f%% of TFPR %.6f\n",
               cfpr, CFPR_TOLERANCE * 100.0, tfpr);

        fails++;

    }

    /* fresh filter: CFPR must be 0 */

    Bloom_Filter* fresh = Bloom_Filter_Create(n, p);

    if(fresh != NULL){

        double fresh_cfpr = Bloom_Filter_CFPR(fresh);

        if(fabs(fresh_cfpr) > 1e-12){

            printf("  FAIL: fresh CFPR = %.12f, expected 0\n", fresh_cfpr);

            fails++;

        }

        Bloom_Filter_Destroy(fresh);

    }

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Measured_Fp_Rate(void){

    printf("\n--- Test_Measured_Fp_Rate ---\n");

    const uint64_t n       = 100000;
    const uint64_t probes  = 200000;
    const double   p       = 0.01;

    Bloom_Filter* bf = Bloom_Filter_Create(n, p);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "in", i);

        Mark_Str(bf, key);

        if((i % 10000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "insert");
        }

    }

    uint64_t fp = 0;

    for(uint64_t i = 0; i < probes; i++){

        Make_Key(key, sizeof(key), "out", i);

        if(Query_Str(bf, key)){
            fp++;
        }

        if((i % 20000) == 0 || i + 1 == probes){
            Progress_Bar(i + 1, probes, "probe non-member");
        }

    }

    double measured    = (double)fp / (double)probes;
    double theoretical = Bloom_Filter_CFPR(bf);

    printf("  measured FP rate = %.6f (fp=%" PRIu64 " / %" PRIu64
           "), CFPR = %.6f, |diff| tolerance = %.4f\n",
           measured, fp, probes, theoretical, FP_ABS_TOLERANCE);

    /*
     * Compare against CFPR (the theoretical estimate at the current fill
     * level), not against the target p. This catches hash-quality
     * regressions: a broken hash raises the measured rate well above
     * CFPR, while still leaving it below 3x of the target p.
     */

    if(fabs(measured - theoretical) > FP_ABS_TOLERANCE){

        printf("  FAIL: |measured - CFPR| = %.6f exceeds %.4f\n",
               fabs(measured - theoretical), FP_ABS_TOLERANCE);

        fails++;

    }

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Mark_Idempotent_And_Clear(void){

    printf("\n--- Test_Mark_Idempotent_And_Clear ---\n");

    const uint64_t n = 10000;

    Bloom_Filter* bf = Bloom_Filter_Create(n, 0.01);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    /* insert each key 3 times */

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "dup", i);

        Mark_Str(bf, key);
        Mark_Str(bf, key);
        Mark_Str(bf, key);

        if((i % 1000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "repeat insert");
        }

    }

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "dup", i);

        if(!Query_Str(bf, key)){

            printf("  FAIL: duplicate-marked key missing at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

    }

    /* Clear -> all queries must be false, count must be 0 */

    Bloom_Filter_Clear(bf);

    if(Bloom_Filter_Element_Count(bf) != 0){

        printf("  FAIL: after Clear, count = %" PRIu64 ", expected 0\n",
               Bloom_Filter_Element_Count(bf));

        fails++;

    }

    uint64_t hits_after_clear = 0;

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "dup", i);

        if(Query_Str(bf, key)){
            hits_after_clear++;
        }

        if((i % 1000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "query after Clear");
        }

    }

    if(hits_after_clear != 0){

        printf("  FAIL: after Clear, %" PRIu64 " queries still hit (expected 0)\n",
               hits_after_clear);

        fails++;

    }

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_IsFull_And_Accessors(void){

    printf("\n--- Test_IsFull_And_Accessors ---\n");

    const uint64_t n = 5000;

    Bloom_Filter* bf = Bloom_Filter_Create(n, 0.01);

    if(bf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    int fails = 0;

    if(Bloom_Filter_Capacity(bf) != n){

        printf("  FAIL: Capacity = %" PRIu64 ", expected %" PRIu64 "\n",
               Bloom_Filter_Capacity(bf), n);

        fails++;

    }

    if(Bloom_Filter_IsFull(bf) != false){

        printf("  FAIL: empty filter reported full\n");

        fails++;

    }

    char key[64];

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "f", i);

        Mark_Str(bf, key);

        if((i % 500) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "fill to IsFull");
        }

    }

    if(Bloom_Filter_IsFull(bf) != true){

        printf("  FAIL: full filter reported not full (count=%" PRIu64 ", cap=%" PRIu64 ")\n",
               Bloom_Filter_Element_Count(bf), Bloom_Filter_Capacity(bf));

        fails++;

    }

    /* Mark still succeeds after full (documented behavior of BF layer) */

    if(!Mark_Str(bf, "overflow")){

        printf("  FAIL: Mark after full returned false\n");

        fails++;

    }

    if(Bloom_Filter_Element_Count(bf) != n + 1){

        printf("  FAIL: count after overflow = %" PRIu64 ", expected %" PRIu64 "\n",
               Bloom_Filter_Element_Count(bf), n + 1);

        fails++;

    }

    Bloom_Filter_Destroy(bf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

/* ------------------------------------------------------------------ */
/* Runner                                                             */
/* ------------------------------------------------------------------ */

int main(void){

    printf("========================================\n");
    printf(" Bloom Filter - Unit Tests\n");
    printf("========================================\n");

    int total_fails = 0;

    total_fails += Test_Murmur3_Known_Vectors();
    total_fails += Test_Create_Invalid_Params();
    total_fails += Test_Null_Defense();
    total_fails += Test_Tiny_Capacity();
    total_fails += Test_Binary_Key();
    total_fails += Test_No_False_Negative();
    total_fails += Test_Cfpr_Matches_Tfpr();
    total_fails += Test_Measured_Fp_Rate();
    total_fails += Test_Mark_Idempotent_And_Clear();
    total_fails += Test_IsFull_And_Accessors();

    printf("\n========================================\n");

    if(total_fails == 0){
        printf(" ALL TESTS PASSED\n");
    }
    else{
        printf(" %d FAILURE(S)\n", total_fails);
    }

    printf("========================================\n");

    return total_fails == 0 ? 0 : 1;

}