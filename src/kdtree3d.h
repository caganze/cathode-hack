/**
 * kdtree3d.h
 * 
 * Fast CPU 3D KD-tree for nearest-neighbor search
 * Adapted from fast-bosquey-tree (https://github.com/caganze/fast-bosquey-tree)
 * 
 * Key features:
 * - In-place tree construction (memory efficient)
 * - O(log N) queries vs O(N) brute force
 * - CPU batch processing with OpenMP parallelization
 * - Bucket size = 32 (cache-line optimal)
 * - Custom scaling factors for anisotropic distance metrics
 * 
 * Usage:
 *   KDTree3D tree;
 *   kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);
 *   kdtree3d_build(x, y, z, n, &tree);
 *   
 *   int nearest_idx = kdtree3d_nearest(&tree, qx, qy, qz, &dist_sq);
 *   
 *   kdtree3d_free(&tree);
 */

#ifndef KDTREE3D_H
#define KDTREE3D_H

#include <stddef.h>

/* Bucket size - leaves contain up to this many points */
/* 32 is optimal for cache lines */
#define KDTREE3D_BUCKET_SIZE 32

/* 3D point with original index */
typedef struct {
    float x, y, z;
    int original_idx;
} KDTree3DPoint;

/* CPU KD-tree structure (flat array representation) */
typedef struct {
    KDTree3DPoint *points;  /* Tree-ordered points */
    int n_points;
    int is_built;
    
    /* Scaling factors for distance calculation */
    float scale_x, scale_y, scale_z;
} KDTree3D;

/* ============================================================================
 * CPU Functions
 * ============================================================================ */

/**
 * Initialize tree with custom scaling factors
 * Scaling is applied during distance calculation: d = (dx*sx)^2 + (dy*sy)^2 + (dz*sz)^2
 * 
 * @param tree: Tree structure to initialize
 * @param scale_x, scale_y, scale_z: Scaling factors (default 1.0)
 */
void kdtree3d_init(KDTree3D *tree, float scale_x, float scale_y, float scale_z);

/**
 * Build KD-tree from 3D points (CPU, in-place)
 * 
 * @param x, y, z: Input coordinate arrays
 * @param n_points: Number of points
 * @param tree: Output tree structure (must be initialized or zeroed)
 * @return 0 on success, -1 on error
 */
int kdtree3d_build(const float *x, const float *y, const float *z, 
                   int n_points, KDTree3D *tree);

/**
 * Query nearest neighbor (CPU, single query)
 * 
 * @param tree: Built tree
 * @param qx, qy, qz: Query point
 * @param dist_sq_out: Output squared distance (optional, can be NULL)
 * @return Index into original array
 */
int kdtree3d_nearest(const KDTree3D *tree, float qx, float qy, float qz, 
                     float *dist_sq_out);

/**
 * Query nearest neighbors (CPU, batch)
 * 
 * @param tree: Built tree
 * @param qx, qy, qz: Query point arrays
 * @param result_idx: Output indices into original array
 * @param n_queries: Number of queries
 */
void kdtree3d_nearest_batch(const KDTree3D *tree,
                            const float *qx, const float *qy, const float *qz,
                            int *result_idx, int n_queries);

/**
 * Query k nearest neighbors (CPU, single query)
 * 
 * @param tree: Built tree
 * @param qx, qy, qz: Query point
 * @param k: Number of neighbors to find
 * @param result_idx: Output array of indices (must be pre-allocated to size k)
 * @param result_dist_sq: Output array of squared distances (optional, can be NULL)
 * @return Actual number of neighbors found (may be < k if tree has fewer points)
 */
int kdtree3d_knearest(const KDTree3D *tree, float qx, float qy, float qz,
                      int k, int *result_idx, float *result_dist_sq);

/**
 * Free tree memory
 */
void kdtree3d_free(KDTree3D *tree);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Get tree statistics
 */
void kdtree3d_stats(const KDTree3D *tree, int *depth_out, int *n_leaves_out);

/**
 * Validate tree structure
 * 
 * @return 0 if valid, -1 if corrupted
 */
int kdtree3d_validate(const KDTree3D *tree);

#endif /* KDTREE3D_H */
