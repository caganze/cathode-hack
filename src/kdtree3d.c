/**
 * kdtree3d.c
 * 
 * Fast CPU 3D KD-tree implementation
 * Adapted from fast-bosquey-tree (https://github.com/caganze/fast-bosquey-tree)
 */

#include "kdtree3d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

static inline void swap_points(KDTree3DPoint *a, KDTree3DPoint *b) {
    KDTree3DPoint tmp = *a;
    *a = *b;
    *b = tmp;
}

static inline float get_dim(const KDTree3DPoint *p, int dim) {
    switch (dim) {
        case 0: return p->x;
        case 1: return p->y;
        case 2: return p->z;
        default: return p->x;
    }
}

static int partition(KDTree3DPoint *points, int left, int right, int pivot_idx, int dim) {
    float pivot_val = get_dim(&points[pivot_idx], dim);
    swap_points(&points[pivot_idx], &points[right]);
    int store_idx = left;
    
    for (int i = left; i < right; i++) {
        if (get_dim(&points[i], dim) < pivot_val) {
            swap_points(&points[store_idx], &points[i]);
            store_idx++;
        }
    }
    swap_points(&points[store_idx], &points[right]);
    return store_idx;
}

static void quickselect(KDTree3DPoint *points, int left, int right, int k, int dim) {
    while (left < right) {
        int pivot_idx = left + (right - left) / 2;
        pivot_idx = partition(points, left, right, pivot_idx, dim);
        
        if (k == pivot_idx) return;
        else if (k < pivot_idx) right = pivot_idx - 1;
        else left = pivot_idx + 1;
    }
}

static void build_recursive(KDTree3DPoint *points, int left, int right, int depth) {
    if (right - left + 1 <= KDTREE3D_BUCKET_SIZE) return;
    
    int dim = depth % 3;
    int median = left + (right - left) / 2;
    
    quickselect(points, left, right, median, dim);
    
    if (median > left) build_recursive(points, left, median - 1, depth + 1);
    if (median < right) build_recursive(points, median + 1, right, depth + 1);
}

/* CPU distance with scaling */
static inline float dist_sq_cpu(float x1, float y1, float z1,
                                float x2, float y2, float z2,
                                float sx, float sy, float sz) {
    float dx = (x1 - x2) * sx;
    float dy = (y1 - y2) * sy;
    float dz = (z1 - z2) * sz;
    return dx*dx + dy*dy + dz*dz;
}

/* ============================================================================
 * Public API Functions
 * ============================================================================ */

void kdtree3d_init(KDTree3D *tree, float scale_x, float scale_y, float scale_z) {
    tree->points = NULL;
    tree->n_points = 0;
    tree->is_built = 0;
    tree->scale_x = scale_x;
    tree->scale_y = scale_y;
    tree->scale_z = scale_z;
}

int kdtree3d_build(const float *x, const float *y, const float *z,
                   int n_points, KDTree3D *tree) {
    /* Set default scaling if not initialized */
    if (tree->scale_x == 0) tree->scale_x = 1.0f;
    if (tree->scale_y == 0) tree->scale_y = 1.0f;
    if (tree->scale_z == 0) tree->scale_z = 1.0f;
    
    size_t mem_size = n_points * sizeof(KDTree3DPoint);
    
    tree->points = (KDTree3DPoint*)malloc(mem_size);
    if (!tree->points) {
        fprintf(stderr, "Error: Cannot allocate KD-tree\n");
        return -1;
    }
    
    /* Copy data with original indices */
    for (int i = 0; i < n_points; i++) {
        tree->points[i].x = x[i];
        tree->points[i].y = y[i];
        tree->points[i].z = z[i];
        tree->points[i].original_idx = i;
    }
    
    /* Build in-place */
    build_recursive(tree->points, 0, n_points - 1, 0);
    
    tree->n_points = n_points;
    tree->is_built = 1;
    
    return 0;
}

int kdtree3d_nearest(const KDTree3D *tree, float qx, float qy, float qz,
                     float *dist_sq_out) {
    if (!tree->is_built) return -1;
    
    float sx = tree->scale_x, sy = tree->scale_y, sz = tree->scale_z;
    float best_dist = FLT_MAX;
    int best_idx = 0;
    
    /* Stack-based traversal */
    int stack_left[64], stack_right[64], stack_depth[64];
    int sp = 0;
    
    stack_left[0] = 0;
    stack_right[0] = tree->n_points - 1;
    stack_depth[0] = 0;
    sp = 1;
    
    while (sp > 0) {
        sp--;
        int left = stack_left[sp];
        int right = stack_right[sp];
        int depth = stack_depth[sp];
        int size = right - left + 1;
        
        /* Leaf: linear search */
        if (size <= KDTREE3D_BUCKET_SIZE) {
            for (int i = left; i <= right; i++) {
                float d = dist_sq_cpu(qx, qy, qz,
                                      tree->points[i].x, tree->points[i].y, tree->points[i].z,
                                      sx, sy, sz);
                if (d < best_dist) {
                    best_dist = d;
                    best_idx = tree->points[i].original_idx;
                }
            }
            continue;
        }
        
        /* Internal node */
        int dim = depth % 3;
        int median = left + (right - left) / 2;
        
        float d = dist_sq_cpu(qx, qy, qz,
                              tree->points[median].x, tree->points[median].y, tree->points[median].z,
                              sx, sy, sz);
        if (d < best_dist) {
            best_dist = d;
            best_idx = tree->points[median].original_idx;
        }
        
        float query_val = (dim == 0) ? qx * sx : (dim == 1) ? qy * sy : qz * sz;
        float split_val = (dim == 0) ? tree->points[median].x * sx : 
                          (dim == 1) ? tree->points[median].y * sy : 
                                       tree->points[median].z * sz;
        float diff = query_val - split_val;
        
        int near_l, near_r, far_l, far_r;
        if (diff <= 0) {
            near_l = left; near_r = median - 1;
            far_l = median + 1; far_r = right;
        } else {
            near_l = median + 1; near_r = right;
            far_l = left; far_r = median - 1;
        }
        
        if (far_l <= far_r && diff * diff < best_dist) {
            stack_left[sp] = far_l;
            stack_right[sp] = far_r;
            stack_depth[sp] = depth + 1;
            sp++;
        }
        if (near_l <= near_r) {
            stack_left[sp] = near_l;
            stack_right[sp] = near_r;
            stack_depth[sp] = depth + 1;
            sp++;
        }
    }
    
    if (dist_sq_out) *dist_sq_out = best_dist;
    return best_idx;
}

void kdtree3d_nearest_batch(const KDTree3D *tree,
                            const float *qx, const float *qy, const float *qz,
                            int *result_idx, int n_queries) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < n_queries; i++) {
        result_idx[i] = kdtree3d_nearest(tree, qx[i], qy[i], qz[i], NULL);
    }
}

/* ============================================================================
 * K-Nearest Neighbors Query
 * ============================================================================ */

/* Simple max-heap for k-nearest neighbors */
typedef struct {
    float dist_sq;
    int idx;
} HeapEntry;

static inline void heap_swap(HeapEntry *a, HeapEntry *b) {
    HeapEntry tmp = *a;
    *a = *b;
    *b = tmp;
}

static void heap_push(HeapEntry *heap, int *heap_size, int k, float dist_sq, int idx) {
    if (*heap_size < k) {
        /* Heap not full, just add */
        heap[*heap_size].dist_sq = dist_sq;
        heap[*heap_size].idx = idx;
        (*heap_size)++;
        
        /* Bubble up */
        int i = *heap_size - 1;
        while (i > 0) {
            int parent = (i - 1) / 2;
            if (heap[i].dist_sq > heap[parent].dist_sq) {
                heap_swap(&heap[i], &heap[parent]);
                i = parent;
            } else {
                break;
            }
        }
    } else if (dist_sq < heap[0].dist_sq) {
        /* Replace root (largest) */
        heap[0].dist_sq = dist_sq;
        heap[0].idx = idx;
        
        /* Bubble down */
        int i = 0;
        while (1) {
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            int largest = i;
            
            if (left < *heap_size && heap[left].dist_sq > heap[largest].dist_sq)
                largest = left;
            if (right < *heap_size && heap[right].dist_sq > heap[largest].dist_sq)
                largest = right;
            
            if (largest != i) {
                heap_swap(&heap[i], &heap[largest]);
                i = largest;
            } else {
                break;
            }
        }
    }
}

int kdtree3d_knearest(const KDTree3D *tree, float qx, float qy, float qz,
                      int k, int *result_idx, float *result_dist_sq) {
    if (!tree->is_built || k <= 0) return 0;
    if (k > tree->n_points) k = tree->n_points;
    
    float sx = tree->scale_x, sy = tree->scale_y, sz = tree->scale_z;
    
    /* Max-heap to track k nearest */
    HeapEntry *heap = (HeapEntry*)malloc(k * sizeof(HeapEntry));
    if (!heap) return 0;
    int heap_size = 0;
    
    /* Stack-based traversal */
    int stack_left[64], stack_right[64], stack_depth[64];
    int sp = 0;
    
    stack_left[0] = 0;
    stack_right[0] = tree->n_points - 1;
    stack_depth[0] = 0;
    sp = 1;
    
    while (sp > 0) {
        sp--;
        int left = stack_left[sp];
        int right = stack_right[sp];
        int depth = stack_depth[sp];
        int size = right - left + 1;
        
        /* Leaf: linear search */
        if (size <= KDTREE3D_BUCKET_SIZE) {
            for (int i = left; i <= right; i++) {
                float d = dist_sq_cpu(qx, qy, qz,
                                      tree->points[i].x, tree->points[i].y, tree->points[i].z,
                                      sx, sy, sz);
                heap_push(heap, &heap_size, k, d, tree->points[i].original_idx);
            }
            continue;
        }
        
        /* Internal node */
        int dim = depth % 3;
        int median = left + (right - left) / 2;
        
        float d = dist_sq_cpu(qx, qy, qz,
                              tree->points[median].x, tree->points[median].y, tree->points[median].z,
                              sx, sy, sz);
        heap_push(heap, &heap_size, k, d, tree->points[median].original_idx);
        
        float query_val = (dim == 0) ? qx * sx : (dim == 1) ? qy * sy : qz * sz;
        float split_val = (dim == 0) ? tree->points[median].x * sx : 
                          (dim == 1) ? tree->points[median].y * sy : 
                                       tree->points[median].z * sz;
        float diff = query_val - split_val;
        
        int near_l, near_r, far_l, far_r;
        if (diff <= 0) {
            near_l = left; near_r = median - 1;
            far_l = median + 1; far_r = right;
        } else {
            near_l = median + 1; near_r = right;
            far_l = left; far_r = median - 1;
        }
        
        /* Only explore far side if it could contain closer points */
        float prune_dist = (heap_size < k) ? FLT_MAX : heap[0].dist_sq;
        if (far_l <= far_r && diff * diff < prune_dist) {
            stack_left[sp] = far_l;
            stack_right[sp] = far_r;
            stack_depth[sp] = depth + 1;
            sp++;
        }
        if (near_l <= near_r) {
            stack_left[sp] = near_l;
            stack_right[sp] = near_r;
            stack_depth[sp] = depth + 1;
            sp++;
        }
    }
    
    /* Extract results (heap is max-heap, so results are in reverse distance order) */
    /* Sort by distance for consistent output */
    for (int i = 0; i < heap_size - 1; i++) {
        for (int j = i + 1; j < heap_size; j++) {
            if (heap[j].dist_sq < heap[i].dist_sq) {
                heap_swap(&heap[i], &heap[j]);
            }
        }
    }
    
    for (int i = 0; i < heap_size; i++) {
        result_idx[i] = heap[i].idx;
        if (result_dist_sq) result_dist_sq[i] = heap[i].dist_sq;
    }
    
    free(heap);
    return heap_size;
}

void kdtree3d_free(KDTree3D *tree) {
    if (tree->points) {
        free(tree->points);
        tree->points = NULL;
    }
    tree->is_built = 0;
}

/* ============================================================================
 * Utilities
 * ============================================================================ */

void kdtree3d_stats(const KDTree3D *tree, int *depth_out, int *n_leaves_out) {
    if (!tree->is_built) {
        if (depth_out) *depth_out = 0;
        if (n_leaves_out) *n_leaves_out = 0;
        return;
    }
    
    int depth = (int)(log2(tree->n_points / KDTREE3D_BUCKET_SIZE) + 1);
    int n_leaves = (tree->n_points + KDTREE3D_BUCKET_SIZE - 1) / KDTREE3D_BUCKET_SIZE;
    
    if (depth_out) *depth_out = depth;
    if (n_leaves_out) *n_leaves_out = n_leaves;
}

int kdtree3d_validate(const KDTree3D *tree) {
    if (!tree->is_built || !tree->points) return -1;
    
    /* Check all original indices are valid */
    for (int i = 0; i < tree->n_points; i++) {
        if (tree->points[i].original_idx < 0 || tree->points[i].original_idx >= tree->n_points) {
            return -1;
        }
    }
    
    return 0;
}
