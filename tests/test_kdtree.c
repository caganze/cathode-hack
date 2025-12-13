/**
 * @file test_kdtree.c
 * @brief Unit tests for the KD-tree implementation
 * 
 * Tests:
 * - Tree construction
 * - Nearest neighbor queries
 * - K-nearest neighbor queries
 * - Edge cases (empty tree, single point, duplicate points)
 * - Scaling factors
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "../src/kdtree3d.h"

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

static float rand_float(float min, float max) {
    return min + (float)rand() / RAND_MAX * (max - min);
}

/* Brute-force nearest neighbor for verification */
static int brute_force_nearest(const float *x, const float *y, const float *z,
                               int n, float qx, float qy, float qz,
                               float sx, float sy, float sz, float *dist_sq) {
    int best_idx = 0;
    float best_dist = 1e30f;
    
    for (int i = 0; i < n; i++) {
        float dx = (x[i] - qx) * sx;
        float dy = (y[i] - qy) * sy;
        float dz = (z[i] - qz) * sz;
        float d = dx*dx + dy*dy + dz*dz;
        if (d < best_dist) {
            best_dist = d;
            best_idx = i;
        }
    }
    
    if (dist_sq) *dist_sq = best_dist;
    return best_idx;
}

/* ============================================================================
 * Test Cases
 * ============================================================================ */

void test_tree_init(void) {
    printf("\n--- Testing tree initialization ---\n");
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    
    test_result("Tree initialized with default scaling", 
                tree.scale_x == 1.0f && tree.scale_y == 1.0f && tree.scale_z == 1.0f);
    test_result("Tree marked as not built", tree.is_built == 0);
    test_result("Tree points is NULL", tree.points == NULL);
    
    kdtree3d_free(&tree);
    
    /* Custom scaling */
    kdtree3d_init(&tree, 2.0f, 3.0f, 0.5f);
    test_result("Custom scaling x", tree.scale_x == 2.0f);
    test_result("Custom scaling y", tree.scale_y == 3.0f);
    test_result("Custom scaling z", tree.scale_z == 0.5f);
    
    kdtree3d_free(&tree);
}

void test_tree_build_simple(void) {
    printf("\n--- Testing simple tree build ---\n");
    
    float x[] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f};
    float y[] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f};
    float z[] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f};
    int n = 5;
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    
    int result = kdtree3d_build(x, y, z, n, &tree);
    
    test_result("Build returns success", result == 0);
    test_result("Tree marked as built", tree.is_built == 1);
    test_result("Tree has correct point count", tree.n_points == n);
    test_result("Tree points allocated", tree.points != NULL);
    
    /* Validate tree structure */
    int valid = kdtree3d_validate(&tree);
    test_result("Tree structure is valid", valid == 0);
    
    kdtree3d_free(&tree);
}

void test_tree_build_empty(void) {
    printf("\n--- Testing empty tree build ---\n");
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    
    int result = kdtree3d_build(NULL, NULL, NULL, 0, &tree);
    
    test_result("Build with 0 points returns error", result == -1);
    test_result("Tree not marked as built", tree.is_built == 0);
    
    kdtree3d_free(&tree);
}

void test_tree_build_single(void) {
    printf("\n--- Testing single-point tree ---\n");
    
    float x[] = {5.0f};
    float y[] = {3.0f};
    float z[] = {7.0f};
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    
    int result = kdtree3d_build(x, y, z, 1, &tree);
    
    test_result("Single-point build succeeds", result == 0);
    test_result("Tree has 1 point", tree.n_points == 1);
    
    /* Query the single point */
    float dist_sq;
    int nearest = kdtree3d_nearest(&tree, 5.0f, 3.0f, 7.0f, &dist_sq);
    
    test_result("Nearest returns index 0", nearest == 0);
    test_result("Distance to self is 0", dist_sq < EPSILON);
    
    kdtree3d_free(&tree);
}

void test_nearest_neighbor_exact(void) {
    printf("\n--- Testing exact nearest neighbor ---\n");
    
    /* Grid of 8 points at corners of unit cube */
    float x[] = {0, 0, 0, 0, 1, 1, 1, 1};
    float y[] = {0, 0, 1, 1, 0, 0, 1, 1};
    float z[] = {0, 1, 0, 1, 0, 1, 0, 1};
    int n = 8;
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree);
    
    /* Query at each corner - should return itself */
    int all_correct = 1;
    for (int i = 0; i < n; i++) {
        float dist_sq;
        int nearest = kdtree3d_nearest(&tree, x[i], y[i], z[i], &dist_sq);
        if (dist_sq > EPSILON) {
            all_correct = 0;
        }
    }
    test_result("All corners find themselves", all_correct);
    
    /* Query at center (0.5, 0.5, 0.5) - all corners equidistant */
    float dist_sq;
    int nearest = kdtree3d_nearest(&tree, 0.5f, 0.5f, 0.5f, &dist_sq);
    test_result("Center query returns valid index", nearest >= 0 && nearest < n);
    test_result("Center distance is correct", fabsf(dist_sq - 0.75f) < EPSILON);
    
    kdtree3d_free(&tree);
}

void test_nearest_neighbor_random(void) {
    printf("\n--- Testing random nearest neighbor ---\n");
    
    srand(12345);
    int n = 1000;
    
    float *x = malloc(n * sizeof(float));
    float *y = malloc(n * sizeof(float));
    float *z = malloc(n * sizeof(float));
    
    for (int i = 0; i < n; i++) {
        x[i] = rand_float(-100.0f, 100.0f);
        y[i] = rand_float(-100.0f, 100.0f);
        z[i] = rand_float(-100.0f, 100.0f);
    }
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree);
    
    /* Test 100 random queries against brute force */
    int matches = 0;
    for (int q = 0; q < 100; q++) {
        float qx = rand_float(-100.0f, 100.0f);
        float qy = rand_float(-100.0f, 100.0f);
        float qz = rand_float(-100.0f, 100.0f);
        
        float kd_dist, bf_dist;
        int kd_idx = kdtree3d_nearest(&tree, qx, qy, qz, &kd_dist);
        int bf_idx = brute_force_nearest(x, y, z, n, qx, qy, qz, 1, 1, 1, &bf_dist);
        
        /* Allow for ties (same distance) */
        if (fabsf(kd_dist - bf_dist) < EPSILON) {
            matches++;
        }
    }
    
    test_result("KD-tree matches brute force for random queries", matches == 100);
    
    free(x); free(y); free(z);
    kdtree3d_free(&tree);
}

void test_knearest(void) {
    printf("\n--- Testing k-nearest neighbors ---\n");
    
    /* Line of points: 0, 1, 2, ..., 9 */
    float x[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    float y[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float z[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int n = 10;
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree);
    
    int result_idx[5];
    float result_dist[5];
    
    /* Query at 5.0 - neighbors should be 5, 4, 6, 3, 7 */
    int found = kdtree3d_knearest(&tree, 5.0f, 0.0f, 0.0f, 5, result_idx, result_dist);
    
    test_result("Found 5 neighbors", found == 5);
    test_result("First neighbor is at distance 0", result_dist[0] < EPSILON);
    
    /* Check that all 5 nearest are in result (order may vary for equidistant) */
    int has_5 = 0, has_4 = 0, has_6 = 0;
    for (int i = 0; i < found; i++) {
        int orig_idx = result_idx[i];
        if (x[orig_idx] == 5.0f) has_5 = 1;
        if (x[orig_idx] == 4.0f) has_4 = 1;
        if (x[orig_idx] == 6.0f) has_6 = 1;
    }
    test_result("k-nearest contains 5", has_5);
    test_result("k-nearest contains 4", has_4);
    test_result("k-nearest contains 6", has_6);
    
    kdtree3d_free(&tree);
}

void test_knearest_more_than_n(void) {
    printf("\n--- Testing k > n neighbors ---\n");
    
    float x[] = {0, 1, 2};
    float y[] = {0, 0, 0};
    float z[] = {0, 0, 0};
    int n = 3;
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree);
    
    int result_idx[10];
    float result_dist[10];
    
    /* Ask for 10 neighbors when only 3 exist */
    int found = kdtree3d_knearest(&tree, 0.0f, 0.0f, 0.0f, 10, result_idx, result_dist);
    
    test_result("Returns min(k, n) neighbors", found == 3);
    
    kdtree3d_free(&tree);
}

void test_scaling_factors(void) {
    printf("\n--- Testing scaling factors ---\n");
    
    /* Two points: (0,0,0) and (1,1,1) */
    float x[] = {0, 1};
    float y[] = {0, 1};
    float z[] = {0, 1};
    int n = 2;
    
    /* With equal scaling, distance^2 = 3 */
    KDTree3D tree1;
    kdtree3d_init(&tree1, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree1);
    
    float dist_sq;
    kdtree3d_nearest(&tree1, 1.0f, 1.0f, 1.0f, &dist_sq);
    test_result("Equal scaling: dist to (1,1,1) is 0", dist_sq < EPSILON);
    
    kdtree3d_nearest(&tree1, 0.0f, 0.0f, 0.0f, &dist_sq);
    test_result("Equal scaling: dist to (0,0,0) is 0", dist_sq < EPSILON);
    kdtree3d_free(&tree1);
    
    /* With z scaled by 0: distance only in x,y plane */
    KDTree3D tree2;
    kdtree3d_init(&tree2, 1.0f, 1.0f, 0.0f);  /* z doesn't count */
    kdtree3d_build(x, y, z, n, &tree2);
    
    /* Query at (0,0,100) - should still find (0,0,0) as nearest */
    int nearest = kdtree3d_nearest(&tree2, 0.0f, 0.0f, 100.0f, &dist_sq);
    test_result("Z-scale 0: finds nearest ignoring z", 
                tree2.points[nearest].x == 0.0f);
    
    kdtree3d_free(&tree2);
}

void test_duplicate_points(void) {
    printf("\n--- Testing duplicate points ---\n");
    
    /* All points at same location */
    float x[] = {5, 5, 5, 5, 5};
    float y[] = {3, 3, 3, 3, 3};
    float z[] = {7, 7, 7, 7, 7};
    int n = 5;
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    
    int result = kdtree3d_build(x, y, z, n, &tree);
    test_result("Build with duplicates succeeds", result == 0);
    
    float dist_sq;
    int nearest = kdtree3d_nearest(&tree, 5.0f, 3.0f, 7.0f, &dist_sq);
    test_result("Query at duplicates finds one", nearest >= 0 && nearest < n);
    test_result("Distance is 0", dist_sq < EPSILON);
    
    kdtree3d_free(&tree);
}

void test_batch_queries(void) {
    printf("\n--- Testing batch queries ---\n");
    
    srand(54321);
    int n = 500;
    
    float *x = malloc(n * sizeof(float));
    float *y = malloc(n * sizeof(float));
    float *z = malloc(n * sizeof(float));
    
    for (int i = 0; i < n; i++) {
        x[i] = rand_float(-10.0f, 10.0f);
        y[i] = rand_float(-10.0f, 10.0f);
        z[i] = rand_float(-10.0f, 10.0f);
    }
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree);
    
    /* Batch of queries */
    int n_queries = 50;
    float *qx = malloc(n_queries * sizeof(float));
    float *qy = malloc(n_queries * sizeof(float));
    float *qz = malloc(n_queries * sizeof(float));
    int *results = malloc(n_queries * sizeof(int));
    
    for (int i = 0; i < n_queries; i++) {
        qx[i] = rand_float(-10.0f, 10.0f);
        qy[i] = rand_float(-10.0f, 10.0f);
        qz[i] = rand_float(-10.0f, 10.0f);
    }
    
    kdtree3d_nearest_batch(&tree, qx, qy, qz, results, n_queries);
    
    /* Verify against single queries */
    int matches = 0;
    for (int i = 0; i < n_queries; i++) {
        float dist_sq;
        int single = kdtree3d_nearest(&tree, qx[i], qy[i], qz[i], &dist_sq);
        if (results[i] == single) matches++;
    }
    
    test_result("Batch matches single queries", matches == n_queries);
    
    free(x); free(y); free(z);
    free(qx); free(qy); free(qz); free(results);
    kdtree3d_free(&tree);
}

void test_tree_stats(void) {
    printf("\n--- Testing tree statistics ---\n");
    
    srand(99999);
    int n = 1000;
    
    float *x = malloc(n * sizeof(float));
    float *y = malloc(n * sizeof(float));
    float *z = malloc(n * sizeof(float));
    
    for (int i = 0; i < n; i++) {
        x[i] = rand_float(0.0f, 100.0f);
        y[i] = rand_float(0.0f, 100.0f);
        z[i] = rand_float(0.0f, 100.0f);
    }
    
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
    kdtree3d_build(x, y, z, n, &tree);
    
    int depth, n_leaves;
    kdtree3d_stats(&tree, &depth, &n_leaves);
    
    test_result("Depth is reasonable", depth > 0 && depth < 20);
    test_result("Has multiple leaves", n_leaves > 1);
    
    printf("    Tree stats: depth=%d, leaves=%d\n", depth, n_leaves);
    
    free(x); free(y); free(z);
    kdtree3d_free(&tree);
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    printf("========================================\n");
    printf("    KD-Tree Unit Tests\n");
    printf("========================================\n");
    
    test_tree_init();
    test_tree_build_simple();
    test_tree_build_empty();
    test_tree_build_single();
    test_nearest_neighbor_exact();
    test_nearest_neighbor_random();
    test_knearest();
    test_knearest_more_than_n();
    test_scaling_factors();
    test_duplicate_points();
    test_batch_queries();
    test_tree_stats();
    
    printf("\n========================================\n");
    printf("    Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");
    
    return (tests_passed == tests_run) ? 0 : 1;
}
