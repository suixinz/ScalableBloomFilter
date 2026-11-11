#include "scalable_bloom_filter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#if defined(_WIN32)
#  include <io.h>
#  define isatty _isatty
#  define fileno _fileno
#else
#  include <unistd.h>
#endif

/* ------------------------------------------------------------------ */
/* Constants                                                          */
/* ------------------------------------------------------------------ */

#define CFPR_SLACK     (1.5)    /* SBF CFPR <= 1.5 * TFPR */
#define FP_ABS_TOL     (0.005)  /* |measured - CFPR| tolerance */

/* ------------------------------------------------------------------ */
/* Test utilities                                                     */
/* ------------------------------------------------------------------ */

static void Progress_Bar(uint64_t step, uint64_t total, const char* label){

    if(total == 0){
        return;
    }

    /*
     * Suppress the bar when stdout is not a terminal (CI logs). The '\r'
     * sequence would otherwise be captured as garbage by the CI log.
     */

    if(!isatty(fileno(stdout))){
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

/* Convenience: mark / query a C string (NUL-terminated). */

static inline bool Mark_Str(SBF* sbf, const char* s){
    return SBF_Element_Mark(sbf, s, (uint64_t)strlen(s));
}

static inline bool Query_Str(const SBF* const sbf, const char* s){
    return SBF_Element_Query(sbf, s, (uint64_t)strlen(s));
}

/* ------------------------------------------------------------------ */
/* Test cases                                                         */
/* ------------------------------------------------------------------ */

static int Test_Create_Invalid(void){

    printf("\n--- Test_Create_Invalid ---\n");

    int fails = 0;

    if(SBF_Create(0.0) != NULL){

        printf("  FAIL: p=0 should fail\n");

        fails++;

    }

    if(SBF_Create(1.0) != NULL){

        printf("  FAIL: p=1 should fail\n");

        fails++;

    }

    if(SBF_Create(-0.1) != NULL){

        printf("  FAIL: p<0 should fail\n");

        fails++;

    }

    if(SBF_Create(1.1) != NULL){

        printf("  FAIL: p>1 should fail\n");

        fails++;

    }

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Null_Defense(void){

    printf("\n--- Test_Null_Defense ---\n");

    int fails = 0;

    SBF_Clear(NULL);
    SBF_Destroy(NULL);

    /* key == NULL with len > 0: invalid. */

    if(SBF_Element_Mark(NULL, "x", 1) != false){

        printf("  FAIL: Mark(NULL, \"x\", 1)\n");

        fails++;

    }

    if(SBF_Element_Query(NULL, "x", 1) != false){

        printf("  FAIL: Query(NULL, \"x\", 1)\n");

        fails++;

    }

    if(SBF_Block_Count(NULL) != UINT64_MAX){

        printf("  FAIL: Block_Count(NULL)\n");

        fails++;

    }

    if(SBF_Ele_Count(NULL) != UINT64_MAX){

        printf("  FAIL: Ele_Count(NULL)\n");

        fails++;

    }

    if(!isnan(SBF_CFPR(NULL))){

        printf("  FAIL: CFPR(NULL)\n");

        fails++;

    }

    if(!isnan(SBF_TFPR(NULL))){

        printf("  FAIL: TFPR(NULL)\n");

        fails++;

    }

    SBF* sbf = SBF_Create(0.01);

    if(sbf == NULL){

        printf("  FAIL: Create failed\n");

        return fails + 1;

    }

    /* key == NULL with len > 0: invalid. */

    if(SBF_Element_Mark(sbf, NULL, 5) != false){

        printf("  FAIL: Mark(sbf, NULL, 5)\n");

        fails++;

    }

    if(SBF_Element_Query(sbf, NULL, 5) != false){

        printf("  FAIL: Query(sbf, NULL, 5)\n");

        fails++;

    }

    /*
     * key == NULL with len == 0 is a valid "empty key": Murmur3 allows
     * (NULL, 0), and the underlying bloom_filter.c accepts it. SBF must
     * inherit that contract.
     */

    if(SBF_Element_Mark(sbf, NULL, 0) != true){

        printf("  FAIL: Mark(sbf, NULL, 0) should succeed (empty key)\n");

        fails++;

    }

    if(SBF_Element_Query(sbf, NULL, 0) != true){

        printf("  FAIL: Query(sbf, NULL, 0) should be true after Mark\n");

        fails++;

    }

    SBF_Destroy(sbf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Cfpr_And_Tfpr(void){

    printf("\n--- Test_Cfpr_And_Tfpr ---\n");

    SBF* sbf = SBF_Create(0.01);

    if(sbf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    int fails = 0;

    if(fabs(SBF_CFPR(sbf)) > 1e-12){

        printf("  FAIL: empty SBF CFPR = %.12f, expected 0\n", SBF_CFPR(sbf));

        fails++;

    }

    if(fabs(SBF_TFPR(sbf) - 0.01) > 1e-12){

        printf("  FAIL: TFPR = %.12f, expected 0.01\n", SBF_TFPR(sbf));

        fails++;

    }

    const uint64_t n = 1000000;

    char key[64];

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "cfpr", i);

        Mark_Str(sbf, key);

        if((i % 100000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "fill first block");
        }

    }

    double cfpr = SBF_CFPR(sbf);
    double tfpr = SBF_TFPR(sbf);

    printf("  TFPR = %.6f, CFPR = %.6f, slack = %.2f\n", tfpr, cfpr, CFPR_SLACK);

    /*
     * CFPR is a *live* estimate: it uses each block's current ele_count.
     * Because actual n_i <= C_i, CFPR <= TFPR. The slack below is a
     * sanity bound against a broken aggregation, not a design target.
     */

    if(cfpr > tfpr * CFPR_SLACK){

        printf("  FAIL: CFPR %.6f exceeds %.2fx TFPR %.6f\n",
               cfpr, CFPR_SLACK, tfpr);

        fails++;

    }

    SBF_Destroy(sbf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_No_False_Negative_Across_Expansion(void){

    printf("\n--- Test_No_False_Negative_Across_Expansion ---\n");

    /* 1.5M > DEFAULT_ELEMENT_CAPACITY (1M) -> forces at least one expansion. */

    const uint64_t n = 1500000;

    SBF* sbf = SBF_Create(0.01);

    if(sbf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    uint64_t blocks_before = SBF_Block_Count(sbf);

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "user", i);

        if(!Mark_Str(sbf, key)){

            printf("  FAIL: Mark failed at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 100000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "insert with expansion");
        }

    }

    uint64_t blocks_after = SBF_Block_Count(sbf);

    printf("  block count: before=%" PRIu64 ", after=%" PRIu64 "\n",
           blocks_before, blocks_after);

    if(blocks_after <= blocks_before){

        printf("  FAIL: expected expansion (block count did not grow)\n");

        fails++;

    }

    for(uint64_t i = 0; i < n && fails == 0; i++){

        Make_Key(key, sizeof(key), "user", i);

        if(!Query_Str(sbf, key)){

            printf("  FAIL: false negative at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 100000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "query across blocks");
        }

    }

    uint64_t ele = SBF_Ele_Count(sbf);

    printf("  Ele_Count = %" PRIu64 " (expected %" PRIu64 ")\n", ele, n);

    if(ele != n){

        printf("  FAIL: element count mismatch\n");

        fails++;

    }

    SBF_Destroy(sbf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Clear_Resets_State(void){

    printf("\n--- Test_Clear_Resets_State ---\n");

    /* Exceed default capacity to force expansion, so Clear must handle >1 block. */

    const uint64_t n = 1200000;

    SBF* sbf = SBF_Create(0.01);

    if(sbf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "pre", i);

        Mark_Str(sbf, key);

        if((i % 100000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "fill before Clear");
        }

    }

    uint64_t blocks_before = SBF_Block_Count(sbf);

    SBF_Clear(sbf);

    if(SBF_Ele_Count(sbf) != 0){

        printf("  FAIL: after Clear, Ele_Count = %" PRIu64 ", expected 0\n",
               SBF_Ele_Count(sbf));

        fails++;

    }

    if(SBF_Block_Count(sbf) != blocks_before){

        printf("  FAIL: Clear must keep block count (got %" PRIu64 ", expected %" PRIu64 ")\n",
               SBF_Block_Count(sbf), blocks_before);

        fails++;

    }

    /* TFPR is a design constant; Clear must not change it. */

    if(fabs(SBF_TFPR(sbf) - 0.01) > 1e-12){

        printf("  FAIL: Clear changed TFPR (got %.12f)\n", SBF_TFPR(sbf));

        fails++;

    }

    uint64_t hits = 0;

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "pre", i);

        if(Query_Str(sbf, key)){
            hits++;
        }

        if((i % 100000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "query after Clear");
        }

    }

    if(hits != 0){

        printf("  FAIL: after Clear, %" PRIu64 " keys still hit (expected 0)\n", hits);

        fails++;

    }

    /* New insertions after Clear must go into head block (cur was reset). */

    if(!Mark_Str(sbf, "post_clear_key")){

        printf("  FAIL: Mark after Clear failed\n");

        fails++;

    }

    if(!Query_Str(sbf, "post_clear_key")){

        printf("  FAIL: Query after Clear missing\n");

        fails++;

    }

    if(SBF_Ele_Count(sbf) != 1){

        printf("  FAIL: after Clear + 1 mark, Ele_Count = %" PRIu64 ", expected 1\n",
               SBF_Ele_Count(sbf));

        fails++;

    }

    SBF_Destroy(sbf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Realistic_User_Registration(void){

    printf("\n--- Test_Realistic_User_Registration ---\n");

    /*
     * Simulate a user-registration dedup filter:
     *   - 800k "registered" users inserted (fits in one block, no expansion)
     *   - query each registered user -> must be true
     *   - query 400k "not registered" users -> measure FP rate
     */

    const uint64_t registered = 800000;
    const uint64_t probes     = 400000;
    const double   p          = 0.01;

    SBF* sbf = SBF_Create(p);

    if(sbf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    for(uint64_t i = 0; i < registered; i++){

        Make_Key(key, sizeof(key), "reg", i);

        if(!Mark_Str(sbf, key)){

            printf("  FAIL: Mark failed at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 50000) == 0 || i + 1 == registered){
            Progress_Bar(i + 1, registered, "register users");
        }

    }

    for(uint64_t i = 0; i < registered && fails == 0; i++){

        Make_Key(key, sizeof(key), "reg", i);

        if(!Query_Str(sbf, key)){

            printf("  FAIL: false negative for registered user i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 50000) == 0 || i + 1 == registered){
            Progress_Bar(i + 1, registered, "verify registered");
        }

    }

    uint64_t fp = 0;

    for(uint64_t i = 0; i < probes; i++){

        Make_Key(key, sizeof(key), "unreg", i);

        if(Query_Str(sbf, key)){
            fp++;
        }

        if((i % 50000) == 0 || i + 1 == probes){
            Progress_Bar(i + 1, probes, "probe unregistered");
        }

    }

    double measured = (double)fp / (double)probes;
    double cfpr     = SBF_CFPR(sbf);

    printf("  measured FP = %.6f (%" PRIu64 "/%" PRIu64 "), CFPR = %.6f, TFPR = %.6f\n",
           measured, fp, probes, cfpr, SBF_TFPR(sbf));

    /*
     * Compare against the live CFPR estimate, not the target P.
     * A broken aggregation (or a broken per-block FPR) raises the
     * measured rate well above CFPR while still leaving it below the
     * 1.5x slack against P. The absolute tolerance catches that.
     */

    if(fabs(measured - cfpr) > FP_ABS_TOL){

        printf("  FAIL: |measured - CFPR| = %.6f exceeds %.4f\n",
               fabs(measured - cfpr), FP_ABS_TOL);

        fails++;

    }

    SBF_Destroy(sbf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

static int Test_Expansion_Stress(void){

    printf("\n--- Test_Expansion_Stress ---\n");

    /* Push through 2+ expansions to exercise the prev-chain query path. */

    const uint64_t n = 3500000;

    SBF* sbf = SBF_Create(0.01);

    if(sbf == NULL){

        printf("  FAIL: Create failed\n");

        return 1;

    }

    char key[64];
    int fails = 0;

    for(uint64_t i = 0; i < n; i++){

        Make_Key(key, sizeof(key), "stress", i);

        if(!Mark_Str(sbf, key)){

            printf("  FAIL: Mark failed at i=%" PRIu64 "\n", i);

            fails++;

            break;

        }

        if((i % 200000) == 0 || i + 1 == n){
            Progress_Bar(i + 1, n, "stress insert");
        }

    }

    uint64_t blocks = SBF_Block_Count(sbf);

    printf("  final block count = %" PRIu64 " (expected >= 3)\n", blocks);

    if(blocks < 3){

        printf("  FAIL: expected at least 3 blocks after 3.5M inserts\n");

        fails++;

    }

    /* Query a sparse sample across all keys to verify prev-chain traversal. */

    const uint64_t samples = 100000;
    uint64_t step = n / samples;

    if(step == 0){
        step = 1;
    }

    const uint64_t total_samples = (n + step - 1) / step;

    for(uint64_t i = 0; i < n && fails == 0; i += step){

        Make_Key(key, sizeof(key), "stress", i);

        if(!Query_Str(sbf, key)){

            printf("  FAIL: false negative at i=%" PRIu64 " (block count=%" PRIu64 ")\n",
                   i, blocks);

            fails++;

            break;

        }

        if((i / step) % 10000 == 0){
            Progress_Bar(i / step + 1, total_samples, "sample query");
        }

    }

    Progress_Bar(total_samples, total_samples, "sample query");

    SBF_Destroy(sbf);

    printf("  %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);

    return fails;

}

/* ------------------------------------------------------------------ */
/* Runner                                                             */
/* ------------------------------------------------------------------ */

int main(void){

    printf("========================================\n");
    printf(" Scalable Bloom Filter - Unit + Integration Tests\n");
    printf("========================================\n");

    int total_fails = 0;

    total_fails += Test_Create_Invalid();
    total_fails += Test_Null_Defense();
    total_fails += Test_Cfpr_And_Tfpr();
    total_fails += Test_No_False_Negative_Across_Expansion();
    total_fails += Test_Clear_Resets_State();
    total_fails += Test_Realistic_User_Registration();
    total_fails += Test_Expansion_Stress();

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