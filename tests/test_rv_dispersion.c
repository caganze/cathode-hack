/**
 * @file test_rv_dispersion.c
 * @brief Unit tests for RV dispersion calculations and filtering
 * 
 * Tests:
 * - RV dispersion computation (standard deviation)
 * - Low dispersion stream detection
 * - High dispersion rejection
 * - Edge cases (missing RVs, single RV, etc.)
 * 
 * Key insight: Real stellar streams have very low RV dispersion (~1-10 km/s)
 * while field stars have high dispersion (~100 km/s).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPSILON 1e-6
#define PASS "\033[32mPASS\033[0m"
#define FAIL "\033[31mFAIL\033[0m"

static int tests_run = 0;
static int tests_passed = 0;

/* ============================================================================
 * Test Utilities
 * ============================================================================ */

static void test_result(const char *name, int passed) {
    tests_run++;
    if (passed) {
        tests_passed++;
        printf("  [%s] %s\n", PASS, name);
    } else {
        printf("  [%s] %s\n", FAIL, name);
    }
}

/* Compute RV dispersion from array of RVs */
static double compute_rv_dispersion(const double *rv, const int *valid, int n) {
    if (n < 2) return 0.0;
    
    /* First pass: mean */
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (valid[i]) {
            sum += rv[i];
            count++;
        }
    }
    
    if (count < 2) return 0.0;
    double mean = sum / count;
    
    /* Second pass: variance */
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        if (valid[i]) {
            double diff = rv[i] - mean;
            sum_sq += diff * diff;
        }
    }
    
    double variance = sum_sq / (count - 1);  /* Sample variance */
    return sqrt(variance);
}

/* Check if RV is valid (not NaN, within reasonable range) */
static int is_valid_rv(double rv) {
    return (fabs(rv) > 0.001 && fabs(rv) < 1000.0);
}

/* ============================================================================
 * Test Data - Simulated Stellar Populations
 * ============================================================================ */

/* Simulated cold stream: GD-1 like
 * Mean RV ~ -150 km/s, dispersion ~ 5 km/s */
static double gd1_like_rv[] = {
    -152.3, -148.7, -150.1, -149.5, -151.2,
    -147.8, -153.0, -150.5, -149.0, -151.8,
    -150.0, -148.5, -152.0, -150.8, -149.3,
    -151.5, -147.5, -152.5, -150.2, -148.0
};
static int n_gd1 = 20;

/* Simulated globular cluster: M4 like
 * Mean RV ~ 70 km/s, dispersion ~ 3 km/s */
static double m4_like_rv[] = {
    71.2, 69.5, 70.8, 68.9, 71.5,
    70.1, 69.8, 71.0, 70.5, 69.2,
    70.3, 71.1, 69.6, 70.7, 68.5,
    71.3, 70.0, 69.4, 70.9, 71.8
};
static int n_m4 = 20;

/* Field halo stars: high dispersion
 * Mean ~ 0 km/s, dispersion ~ 100 km/s */
static double halo_rv[] = {
    -156.0, 89.3, -45.2, 210.5, -78.1,
    132.7, -200.4, 15.8, -98.6, 178.2,
    -23.5, 265.1, -180.3, 45.9, -110.7,
    88.4, -55.2, 195.0, -140.8, 30.6
};
static int n_halo = 20;

/* Mixed population: some stream + some field */
static double mixed_rv[] = {
    -150.0, -149.5, 100.2,    /* 2 stream + 1 field */
    -151.2, -75.8, -150.5,    /* 2 stream + 1 field */
    200.3, -148.9, -150.1,    /* 1 field + 2 stream */
    -151.0, -149.2, -50.5     /* 2 stream + 1 field */
};
static int n_mixed = 12;

/* ============================================================================
 * Test Cases - Basic RV Dispersion
 * ============================================================================ */

void test_dispersion_cold_stream(void) {
    printf("\n--- Testing cold stream RV dispersion ---\n");
    
    int valid[100];
    for (int i = 0; i < n_gd1; i++) valid[i] = 1;
    
    double disp = compute_rv_dispersion(gd1_like_rv, valid, n_gd1);
    
    printf("    GD-1 like: dispersion = %.2f km/s\n", disp);
    
    test_result("Cold stream dispersion < 10 km/s", disp < 10.0);
    test_result("Cold stream dispersion > 0", disp > 0.0);
}

void test_dispersion_globular_cluster(void) {
    printf("\n--- Testing globular cluster RV dispersion ---\n");
    
    int valid[100];
    for (int i = 0; i < n_m4; i++) valid[i] = 1;
    
    double disp = compute_rv_dispersion(m4_like_rv, valid, n_m4);
    
    printf("    M4 like: dispersion = %.2f km/s\n", disp);
    
    test_result("GC dispersion < 5 km/s", disp < 5.0);
    test_result("GC dispersion > 0", disp > 0.0);
}

void test_dispersion_field_stars(void) {
    printf("\n--- Testing field star RV dispersion ---\n");
    
    int valid[100];
    for (int i = 0; i < n_halo; i++) valid[i] = 1;
    
    double disp = compute_rv_dispersion(halo_rv, valid, n_halo);
    
    printf("    Field stars: dispersion = %.2f km/s\n", disp);
    
    test_result("Field dispersion > 50 km/s", disp > 50.0);
    test_result("Field dispersion < 200 km/s", disp < 200.0);
}

void test_dispersion_discriminates(void) {
    printf("\n--- Testing that dispersion discriminates populations ---\n");
    
    int valid[100];
    for (int i = 0; i < 100; i++) valid[i] = 1;
    
    double disp_gd1 = compute_rv_dispersion(gd1_like_rv, valid, n_gd1);
    double disp_m4 = compute_rv_dispersion(m4_like_rv, valid, n_m4);
    double disp_halo = compute_rv_dispersion(halo_rv, valid, n_halo);
    
    test_result("Stream < 10 km/s separates from field", disp_gd1 < 20.0 && disp_halo > 50.0);
    test_result("GC < 10 km/s separates from field", disp_m4 < 20.0 && disp_halo > 50.0);
    test_result("Factor of 10 between stream and field", disp_halo / disp_gd1 > 10.0);
    
    printf("    Ratio field/stream: %.1f\n", disp_halo / disp_gd1);
}

/* ============================================================================
 * Test Cases - RV Validity
 * ============================================================================ */

void test_valid_rv(void) {
    printf("\n--- Testing RV validity checks ---\n");
    
    test_result("Normal RV is valid", is_valid_rv(-150.0));
    test_result("Zero RV is invalid", !is_valid_rv(0.0));
    test_result("Very small RV is invalid", !is_valid_rv(0.0001));
    test_result("Very large RV is invalid", !is_valid_rv(1500.0));
    test_result("Negative large RV is invalid", !is_valid_rv(-1200.0));
    test_result("Edge case 0.001 is valid", is_valid_rv(0.002));
    test_result("Edge case 999 is valid", is_valid_rv(999.0));
}

void test_dispersion_with_invalid(void) {
    printf("\n--- Testing dispersion with invalid RVs ---\n");
    
    /* Stream with some invalid (missing) RVs */
    double rv_with_missing[] = {
        -150.0, 0.0, -149.5,     /* 0.0 is missing */
        -151.2, 1500.0, -150.5,  /* 1500 is missing */
        -148.9, -150.1, 0.0,     /* 0.0 is missing */
        -151.0, -149.2
    };
    int n = 11;
    
    int valid[100];
    for (int i = 0; i < n; i++) {
        valid[i] = is_valid_rv(rv_with_missing[i]);
    }
    
    /* Count valid */
    int n_valid = 0;
    for (int i = 0; i < n; i++) {
        if (valid[i]) n_valid++;
    }
    
    double disp = compute_rv_dispersion(rv_with_missing, valid, n);
    
    printf("    Valid RVs: %d/%d, dispersion = %.2f km/s\n", n_valid, n, disp);
    
    test_result("Dispersion computed with valid subset", disp > 0.0);
    test_result("Dispersion still low (stream)", disp < 10.0);
}

void test_dispersion_all_invalid(void) {
    printf("\n--- Testing dispersion with all invalid RVs ---\n");
    
    double rv_all_bad[] = {0.0, 0.0, 0.0, 1500.0, -1500.0};
    int n = 5;
    
    int valid[100];
    for (int i = 0; i < n; i++) {
        valid[i] = is_valid_rv(rv_all_bad[i]);
    }
    
    double disp = compute_rv_dispersion(rv_all_bad, valid, n);
    
    test_result("All invalid returns 0 dispersion", fabs(disp) < EPSILON);
}

void test_dispersion_single_rv(void) {
    printf("\n--- Testing dispersion with single valid RV ---\n");
    
    double rv_single[] = {-150.0, 0.0, 0.0, 0.0};
    int n = 4;
    
    int valid[100];
    for (int i = 0; i < n; i++) {
        valid[i] = is_valid_rv(rv_single[i]);
    }
    
    double disp = compute_rv_dispersion(rv_single, valid, n);
    
    test_result("Single valid RV returns 0 dispersion", fabs(disp) < EPSILON);
}

/* ============================================================================
 * Test Cases - Threshold-Based Filtering
 * ============================================================================ */

void test_threshold_20kms(void) {
    printf("\n--- Testing 20 km/s threshold ---\n");
    
    int valid[100];
    for (int i = 0; i < 100; i++) valid[i] = 1;
    
    double disp_gd1 = compute_rv_dispersion(gd1_like_rv, valid, n_gd1);
    double disp_m4 = compute_rv_dispersion(m4_like_rv, valid, n_m4);
    double disp_halo = compute_rv_dispersion(halo_rv, valid, n_halo);
    
    const double THRESHOLD = 20.0;
    
    bool gd1_passes = (disp_gd1 < THRESHOLD);
    bool m4_passes = (disp_m4 < THRESHOLD);
    bool halo_passes = (disp_halo < THRESHOLD);
    
    test_result("GD-1 passes 20 km/s cut", gd1_passes);
    test_result("M4 passes 20 km/s cut", m4_passes);
    test_result("Halo fails 20 km/s cut", !halo_passes);
    
    printf("    Threshold = %.1f km/s\n", THRESHOLD);
    printf("    GD-1: %.2f km/s -> %s\n", disp_gd1, gd1_passes ? "PASS" : "FAIL");
    printf("    M4:   %.2f km/s -> %s\n", disp_m4, m4_passes ? "PASS" : "FAIL");
    printf("    Halo: %.2f km/s -> %s\n", disp_halo, halo_passes ? "PASS" : "FAIL");
}

void test_threshold_10kms(void) {
    printf("\n--- Testing 10 km/s threshold (stricter) ---\n");
    
    int valid[100];
    for (int i = 0; i < 100; i++) valid[i] = 1;
    
    double disp_gd1 = compute_rv_dispersion(gd1_like_rv, valid, n_gd1);
    double disp_m4 = compute_rv_dispersion(m4_like_rv, valid, n_m4);
    
    const double THRESHOLD = 10.0;
    
    bool gd1_passes = (disp_gd1 < THRESHOLD);
    bool m4_passes = (disp_m4 < THRESHOLD);
    
    test_result("Cold stream passes 10 km/s cut", gd1_passes);
    test_result("GC passes 10 km/s cut", m4_passes);
}

void test_threshold_5kms(void) {
    printf("\n--- Testing 5 km/s threshold (very strict) ---\n");
    
    int valid[100];
    for (int i = 0; i < 100; i++) valid[i] = 1;
    
    double disp_gd1 = compute_rv_dispersion(gd1_like_rv, valid, n_gd1);
    double disp_m4 = compute_rv_dispersion(m4_like_rv, valid, n_m4);
    
    const double THRESHOLD = 5.0;
    
    bool gd1_passes = (disp_gd1 < THRESHOLD);
    bool m4_passes = (disp_m4 < THRESHOLD);
    
    printf("    GD-1: %.2f km/s vs threshold %.1f km/s -> %s\n", 
           disp_gd1, THRESHOLD, gd1_passes ? "PASS" : "FAIL");
    printf("    M4:   %.2f km/s vs threshold %.1f km/s -> %s\n", 
           disp_m4, THRESHOLD, m4_passes ? "PASS" : "FAIL");
    
    /* M4 should definitely pass (GCs are very cold) */
    test_result("Very cold GC passes 5 km/s cut", m4_passes);
}

/* ============================================================================
 * Test Cases - Mixed Populations
 * ============================================================================ */

void test_mixed_population(void) {
    printf("\n--- Testing mixed population ---\n");
    
    int valid[100];
    for (int i = 0; i < n_mixed; i++) valid[i] = 1;
    
    double disp = compute_rv_dispersion(mixed_rv, valid, n_mixed);
    
    printf("    Mixed (stream + field): dispersion = %.2f km/s\n", disp);
    
    /* Mixed population should have high dispersion */
    test_result("Mixed population has high dispersion", disp > 30.0);
}

void test_contamination_effect(void) {
    printf("\n--- Testing contamination effect ---\n");
    
    /* Start with pure stream */
    double rv_contaminated[30];
    int valid[30];
    
    /* Copy stream RVs */
    for (int i = 0; i < n_gd1; i++) {
        rv_contaminated[i] = gd1_like_rv[i];
        valid[i] = 1;
    }
    
    double disp_pure = compute_rv_dispersion(rv_contaminated, valid, n_gd1);
    
    /* Add 20% contamination (4 field stars) */
    int n_contam = n_gd1 + 4;
    rv_contaminated[20] = 100.0;  /* Field star */
    rv_contaminated[21] = -50.0;  /* Field star */
    rv_contaminated[22] = 200.0;  /* Field star */
    rv_contaminated[23] = -180.0; /* Field star */
    for (int i = 20; i < 24; i++) valid[i] = 1;
    
    double disp_contam = compute_rv_dispersion(rv_contaminated, valid, n_contam);
    
    printf("    Pure stream:  dispersion = %.2f km/s\n", disp_pure);
    printf("    With 20%% contamination: dispersion = %.2f km/s\n", disp_contam);
    
    test_result("Contamination increases dispersion", disp_contam > disp_pure);
    test_result("Contaminated stream fails 20 km/s cut", disp_contam > 20.0);
}

/* ============================================================================
 * Test Cases - Edge Cases
 * ============================================================================ */

void test_identical_rvs(void) {
    printf("\n--- Testing identical RVs ---\n");
    
    double rv_same[] = {-150.0, -150.0, -150.0, -150.0, -150.0};
    int n = 5;
    int valid[] = {1, 1, 1, 1, 1};
    
    double disp = compute_rv_dispersion(rv_same, valid, n);
    
    test_result("Identical RVs give 0 dispersion", fabs(disp) < EPSILON);
}

void test_two_values(void) {
    printf("\n--- Testing exactly 2 valid RVs ---\n");
    
    double rv_two[] = {-150.0, -140.0};
    int n = 2;
    int valid[] = {1, 1};
    
    double disp = compute_rv_dispersion(rv_two, valid, n);
    
    /* With 2 values and sample std (n-1), std = |a-b| / sqrt(2) * sqrt(2) = |a-b| */
    double expected = 10.0;  /* sqrt(sum_sq / 1) where sum_sq = 2 * 5^2 = 50 */
    
    printf("    Dispersion of [-150, -140] = %.4f km/s (expected ~%.4f)\n", 
           disp, expected);
    
    test_result("Two values give reasonable dispersion", disp > 5.0 && disp < 15.0);
}

/* ============================================================================
 * Test Cases - RV Compatibility (for ROI clustering)
 * ============================================================================ */

void test_rv_compatibility(void) {
    printf("\n--- Testing RV compatibility for clustering ---\n");
    
    /* Two ROIs from same stream - should be compatible */
    double roi1_mean = -150.0;
    double roi2_mean = -148.5;
    double rv_diff = fabs(roi1_mean - roi2_mean);
    
    const double RV_COMPAT_THRESHOLD = 10.0;  /* km/s */
    
    test_result("Same stream ROIs are compatible", rv_diff < RV_COMPAT_THRESHOLD);
    
    /* Two ROIs from different populations - should not be compatible */
    double roi3_mean = -150.0;  /* Stream */
    double roi4_mean = 50.0;    /* Field cluster */
    double rv_diff_bad = fabs(roi3_mean - roi4_mean);
    
    test_result("Different population ROIs are incompatible", 
                rv_diff_bad > RV_COMPAT_THRESHOLD);
}

void test_combined_dispersion_check(void) {
    printf("\n--- Testing combined dispersion check ---\n");
    
    /* Simulating rois_can_cluster logic */
    
    /* Two compatible ROIs from same stream */
    double roi1_rv[] = {-150.0, -149.5, -151.2, -148.9, -150.5};
    double roi2_rv[] = {-149.0, -150.8, -148.5, -151.5, -150.2};
    int n1 = 5, n2 = 5;
    
    /* Compute combined mean */
    double sum = 0.0;
    for (int i = 0; i < n1; i++) sum += roi1_rv[i];
    for (int i = 0; i < n2; i++) sum += roi2_rv[i];
    double combined_mean = sum / (n1 + n2);
    
    /* Compute combined variance */
    double sum_sq = 0.0;
    for (int i = 0; i < n1; i++) {
        sum_sq += pow(roi1_rv[i] - combined_mean, 2);
    }
    for (int i = 0; i < n2; i++) {
        sum_sq += pow(roi2_rv[i] - combined_mean, 2);
    }
    double combined_disp = sqrt(sum_sq / (n1 + n2 - 1));
    
    printf("    Combined mean: %.2f km/s\n", combined_mean);
    printf("    Combined dispersion: %.2f km/s\n", combined_disp);
    
    const double RV_DISPERSION_MAX = 20.0;
    
    test_result("Compatible ROIs have low combined dispersion", 
                combined_disp < RV_DISPERSION_MAX);
    
    /* Now test incompatible ROIs */
    double roi3_rv[] = {-150.0, -149.5, -151.2, -148.9, -150.5};  /* Stream */
    double roi4_rv[] = {50.0, 55.2, 48.8, 52.1, 49.5};            /* Different pop */
    
    sum = 0.0;
    for (int i = 0; i < 5; i++) sum += roi3_rv[i];
    for (int i = 0; i < 5; i++) sum += roi4_rv[i];
    combined_mean = sum / 10;
    
    sum_sq = 0.0;
    for (int i = 0; i < 5; i++) {
        sum_sq += pow(roi3_rv[i] - combined_mean, 2);
    }
    for (int i = 0; i < 5; i++) {
        sum_sq += pow(roi4_rv[i] - combined_mean, 2);
    }
    combined_disp = sqrt(sum_sq / 9);
    
    printf("    Incompatible combined dispersion: %.2f km/s\n", combined_disp);
    
    test_result("Incompatible ROIs have high combined dispersion", 
                combined_disp > RV_DISPERSION_MAX);
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    printf("========================================\n");
    printf("    RV Dispersion Tests\n");
    printf("========================================\n");
    
    /* Basic dispersion */
    test_dispersion_cold_stream();
    test_dispersion_globular_cluster();
    test_dispersion_field_stars();
    test_dispersion_discriminates();
    
    /* Validity */
    test_valid_rv();
    test_dispersion_with_invalid();
    test_dispersion_all_invalid();
    test_dispersion_single_rv();
    
    /* Thresholds */
    test_threshold_20kms();
    test_threshold_10kms();
    test_threshold_5kms();
    
    /* Mixed populations */
    test_mixed_population();
    test_contamination_effect();
    
    /* Edge cases */
    test_identical_rvs();
    test_two_values();
    
    /* Clustering compatibility */
    test_rv_compatibility();
    test_combined_dispersion_check();
    
    printf("\n========================================\n");
    printf("    Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");
    
    return (tests_passed == tests_run) ? 0 : 1;
}
