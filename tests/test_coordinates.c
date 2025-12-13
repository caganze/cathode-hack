/**
 * @file test_coordinates.c
 * @brief Unit tests for coordinate transformation functions
 * 
 * Tests:
 * - Equatorial <-> Galactic coordinate transformations
 * - Patch-local coordinate transformations
 * - Proper motion transformations
 * - Angular distance calculations
 * - Edge cases (poles, wrap-around)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Include the main header */
#include "../src/stream_detect.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPSILON 1e-6
#define ANGLE_EPSILON 1e-4  /* ~0.3 arcsec */
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

static double angle_diff(double a1, double a2) {
    /* Compute smallest angle difference, handling wrap-around */
    double diff = fabs(a1 - a2);
    if (diff > 180.0) diff = 360.0 - diff;
    return diff;
}

/* ============================================================================
 * Known Reference Points for Testing
 * ============================================================================ */

/* Galactic North Pole in equatorial coordinates */
#define GNP_RA  192.85948
#define GNP_DEC 27.12825

/* Galactic Center in equatorial coordinates */
#define GC_RA   266.405
#define GC_DEC  -28.936

/* ============================================================================
 * Test Cases - Equatorial/Galactic Transformations
 * ============================================================================ */

void test_galactic_north_pole(void) {
    printf("\n--- Testing Galactic North Pole ---\n");
    
    double l, b;
    equatorial_to_galactic(GNP_RA, GNP_DEC, &l, &b);
    
    /* At the North Galactic Pole, b should be +90 */
    test_result("NGP b = +90 deg", fabs(b - 90.0) < ANGLE_EPSILON);
    
    printf("    NGP: (RA=%.4f, Dec=%.4f) -> (l=%.4f, b=%.4f)\n", 
           GNP_RA, GNP_DEC, l, b);
}

void test_galactic_center(void) {
    printf("\n--- Testing Galactic Center ---\n");
    
    double l, b;
    equatorial_to_galactic(GC_RA, GC_DEC, &l, &b);
    
    /* At the Galactic Center, (l, b) should be ~(0, 0) */
    test_result("GC l ~ 0 deg", fabs(l) < 1.0 || fabs(l - 360.0) < 1.0);
    test_result("GC b ~ 0 deg", fabs(b) < 1.0);
    
    printf("    GC: (RA=%.4f, Dec=%.4f) -> (l=%.4f, b=%.4f)\n", 
           GC_RA, GC_DEC, l, b);
}

void test_round_trip_eq_gal(void) {
    printf("\n--- Testing equatorial <-> galactic round-trip ---\n");
    
    /* Test multiple points */
    double test_points[][2] = {
        {0.0, 0.0},
        {180.0, 45.0},
        {270.0, -30.0},
        {45.0, 89.0},      /* Near north pole */
        {200.0, -89.0},    /* Near south pole */
        {359.999, 0.0},    /* Near RA wrap */
    };
    int n_tests = sizeof(test_points) / sizeof(test_points[0]);
    
    int all_pass = 1;
    for (int i = 0; i < n_tests; i++) {
        double ra_in = test_points[i][0];
        double dec_in = test_points[i][1];
        
        double l, b;
        equatorial_to_galactic(ra_in, dec_in, &l, &b);
        
        double ra_out, dec_out;
        galactic_to_equatorial(l, b, &ra_out, &dec_out);
        
        double ra_diff = angle_diff(ra_in, ra_out);
        double dec_diff = fabs(dec_in - dec_out);
        
        if (ra_diff > ANGLE_EPSILON || dec_diff > ANGLE_EPSILON) {
            all_pass = 0;
            printf("    FAIL: (%.4f, %.4f) -> (%.4f, %.4f) -> (%.4f, %.4f)\n",
                   ra_in, dec_in, l, b, ra_out, dec_out);
        }
    }
    
    test_result("Round-trip preserves coordinates", all_pass);
}

void test_galactic_plane(void) {
    printf("\n--- Testing Galactic Plane points ---\n");
    
    /* Points along the galactic plane should have b ~ 0 */
    double plane_points[][2] = {
        {266.405, -28.936},   /* Galactic center */
        {86.405, 28.936},     /* Galactic anti-center */
    };
    
    int all_pass = 1;
    for (int i = 0; i < 2; i++) {
        double l, b;
        equatorial_to_galactic(plane_points[i][0], plane_points[i][1], &l, &b);
        if (fabs(b) > 2.0) {  /* Allow some tolerance */
            all_pass = 0;
        }
    }
    
    test_result("Galactic plane has b ~ 0", all_pass);
}

/* ============================================================================
 * Test Cases - Angular Distance
 * ============================================================================ */

void test_angular_separation_same_point(void) {
    printf("\n--- Testing angular separation (same point) ---\n");
    
    double sep = angular_separation(180.0, 45.0, 180.0, 45.0);
    test_result("Same point has 0 separation", fabs(sep) < EPSILON);
}

void test_angular_separation_poles(void) {
    printf("\n--- Testing angular separation (poles) ---\n");
    
    /* North to south pole = 180 degrees */
    double sep = angular_separation(0.0, 90.0, 0.0, -90.0);
    test_result("N to S pole = 180 deg", fabs(sep - 180.0) < ANGLE_EPSILON);
    
    /* Equator to north pole = 90 degrees */
    sep = angular_separation(0.0, 0.0, 0.0, 90.0);
    test_result("Equator to N pole = 90 deg", fabs(sep - 90.0) < ANGLE_EPSILON);
}

void test_angular_separation_wrap(void) {
    printf("\n--- Testing angular separation (RA wrap) ---\n");
    
    /* Points at RA=359 and RA=1 should be close (2 degrees apart on equator) */
    double sep = angular_separation(359.0, 0.0, 1.0, 0.0);
    test_result("RA wrap: 359 to 1 = 2 deg", fabs(sep - 2.0) < ANGLE_EPSILON);
    
    /* Test at different declinations */
    sep = angular_separation(359.0, 45.0, 1.0, 45.0);
    double expected = 2.0 * cos(45.0 * M_PI / 180.0);  /* Smaller due to dec */
    test_result("RA wrap at dec=45", fabs(sep - expected) < 0.1);
}

void test_angular_separation_known(void) {
    printf("\n--- Testing angular separation (known distances) ---\n");
    
    /* 90 degrees apart on equator */
    double sep = angular_separation(0.0, 0.0, 90.0, 0.0);
    test_result("90 deg along equator", fabs(sep - 90.0) < ANGLE_EPSILON);
    
    /* Small separation */
    sep = angular_separation(100.0, 50.0, 100.1, 50.0);
    double expected = 0.1 * cos(50.0 * M_PI / 180.0);
    test_result("Small separation at dec=50", fabs(sep - expected) < 0.01);
}

/* ============================================================================
 * Test Cases - Patch Local Coordinates
 * ============================================================================ */

void test_patch_center_transform(void) {
    printf("\n--- Testing patch center transformation ---\n");
    
    Star star;
    star.ra = 180.0;
    star.dec = 45.0;
    
    transform_to_patch_coords(&star, 180.0, 45.0);
    
    /* Star at patch center should have phi=0, lambda=0 */
    test_result("Patch center phi ~ 0", fabs(star.phi) < ANGLE_EPSILON);
    test_result("Patch center lambda ~ 0", fabs(star.lambda) < ANGLE_EPSILON);
}

void test_patch_offset_transform(void) {
    printf("\n--- Testing patch offset transformation ---\n");
    
    Star star;
    star.ra = 181.0;  /* 1 degree offset in RA */
    star.dec = 45.0;
    
    transform_to_patch_coords(&star, 180.0, 45.0);
    
    /* phi should be ~1 degree (corrected for declination) */
    double expected_phi = 1.0 * cos(45.0 * M_PI / 180.0);
    test_result("RA offset gives phi offset", fabs(star.phi - expected_phi) < 0.1);
    test_result("Same dec gives lambda ~ 0", fabs(star.lambda) < ANGLE_EPSILON);
}

void test_patch_dec_offset_transform(void) {
    printf("\n--- Testing patch declination offset ---\n");
    
    Star star;
    star.ra = 180.0;
    star.dec = 46.0;  /* 1 degree offset in Dec */
    
    transform_to_patch_coords(&star, 180.0, 45.0);
    
    /* lambda should be ~1 degree */
    test_result("Dec offset gives lambda ~ 1", fabs(star.lambda - 1.0) < 0.1);
    test_result("Same RA gives phi ~ 0", fabs(star.phi) < 0.1);
}

/* ============================================================================
 * Test Cases - Proper Motion Transformation
 * ============================================================================ */

void test_proper_motion_transform_center(void) {
    printf("\n--- Testing PM transformation at center ---\n");
    
    Star star;
    star.ra = 180.0;
    star.dec = 45.0;
    star.pmra = 5.0;
    star.pmdec = -3.0;
    
    transform_proper_motions(&star, 180.0, 45.0);
    
    /* At patch center, PMs should be similar to input */
    test_result("mu_phi similar to pmra", fabs(star.mu_phi - star.pmra) < 1.0);
    test_result("mu_lambda similar to pmdec", fabs(star.mu_lambda - star.pmdec) < 1.0);
}

void test_proper_motion_magnitude_preserved(void) {
    printf("\n--- Testing PM magnitude preservation ---\n");
    
    Star star;
    star.ra = 190.0;
    star.dec = 40.0;
    star.pmra = 10.0;
    star.pmdec = -5.0;
    
    double pm_total_before = sqrt(star.pmra * star.pmra + star.pmdec * star.pmdec);
    
    transform_proper_motions(&star, 180.0, 45.0);
    
    double pm_total_after = sqrt(star.mu_phi * star.mu_phi + 
                                  star.mu_lambda * star.mu_lambda);
    
    /* Total PM magnitude should be preserved (approximately) */
    double ratio = pm_total_after / pm_total_before;
    test_result("PM magnitude preserved (within 20%)", 
                ratio > 0.8 && ratio < 1.2);
    
    printf("    PM total: %.4f -> %.4f (ratio: %.4f)\n", 
           pm_total_before, pm_total_after, ratio);
}

/* ============================================================================
 * Test Cases - Statistical Functions
 * ============================================================================ */

void test_mean(void) {
    printf("\n--- Testing mean function ---\n");
    
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double m = mean(data, 5);
    
    test_result("Mean of [1,2,3,4,5] = 3", fabs(m - 3.0) < EPSILON);
    
    double data2[] = {10.0, 20.0};
    m = mean(data2, 2);
    test_result("Mean of [10,20] = 15", fabs(m - 15.0) < EPSILON);
}

void test_variance(void) {
    printf("\n--- Testing variance function ---\n");
    
    double data[] = {2.0, 4.0, 6.0, 8.0, 10.0};
    double v = variance(data, 5);
    
    /* Variance of [2,4,6,8,10] = 8 (population variance) */
    test_result("Variance computed correctly", fabs(v - 8.0) < EPSILON);
}

void test_std_dev(void) {
    printf("\n--- Testing standard deviation ---\n");
    
    double data[] = {2.0, 4.0, 6.0, 8.0, 10.0};
    double s = std_dev(data, 5);
    
    test_result("Std dev is sqrt(variance)", fabs(s - sqrt(8.0)) < EPSILON);
}

void test_median(void) {
    printf("\n--- Testing median function ---\n");
    
    double data1[] = {1.0, 3.0, 5.0, 7.0, 9.0};
    double m = median(data1, 5);
    test_result("Median of odd count", fabs(m - 5.0) < EPSILON);
    
    double data2[] = {1.0, 2.0, 3.0, 4.0};
    m = median(data2, 4);
    test_result("Median of even count", fabs(m - 2.5) < EPSILON);
    
    /* Unsorted data */
    double data3[] = {9.0, 1.0, 5.0, 3.0, 7.0};
    m = median(data3, 5);
    test_result("Median of unsorted data", fabs(m - 5.0) < EPSILON);
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    printf("========================================\n");
    printf("    Coordinate Transformation Tests\n");
    printf("========================================\n");
    
    /* Equatorial/Galactic */
    test_galactic_north_pole();
    test_galactic_center();
    test_round_trip_eq_gal();
    test_galactic_plane();
    
    /* Angular separation */
    test_angular_separation_same_point();
    test_angular_separation_poles();
    test_angular_separation_wrap();
    test_angular_separation_known();
    
    /* Patch coordinates */
    test_patch_center_transform();
    test_patch_offset_transform();
    test_patch_dec_offset_transform();
    
    /* Proper motions */
    test_proper_motion_transform_center();
    test_proper_motion_magnitude_preserved();
    
    /* Statistics */
    test_mean();
    test_variance();
    test_std_dev();
    test_median();
    
    printf("\n========================================\n");
    printf("    Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");
    
    return (tests_passed == tests_run) ? 0 : 1;
}
