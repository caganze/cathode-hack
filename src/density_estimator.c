/**
 * @file density_estimator.c
 * @brief Density estimation methods for anomaly detection
 * 
 * Implements various density estimation techniques used to compute
 * anomaly scores R(x) = p_data(x) / p_background(x)
 * 
 * The paper uses normalizing flows, but we provide simpler alternatives:
 * - Kernel Density Estimation (KDE)
 * - K-Nearest Neighbors (KNN) density
 * - Histogram-based density
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "stream_detect.h"
#include "kdtree3d.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Kernel Density Estimation
 * ============================================================================ */

/* 
 * Minimum bandwidths for different dimensions to smooth over survey artifacts.
 * 
 * DESI tiles are ~1.6 deg radius with ~1.4 deg spacing. To properly smooth
 * over the hexagonal tile pattern, we need bandwidths > tile spacing.
 * 
 * - Sky position (phi, lambda): min 2.0 deg to smooth over tile artifacts
 * - Proper motion: no minimum (use Silverman's rule)
 * - Distance: no minimum
 * - RV: no minimum
 */
#define MIN_BANDWIDTH_SKY_DEG    2.0    /* Minimum for phi, lambda dimensions */
#define MIN_BANDWIDTH_PM_MASYR   0.0    /* No minimum for PM (smoother in data) */

/**
 * @brief Create a new KDE estimator with per-dimension bandwidths
 * 
 * @param n_dims Number of dimensions
 * @return Newly allocated KDE, or NULL on error
 */
KDE *kde_create(uint32_t n_dims) {
    KDE *kde = (KDE *)calloc(1, sizeof(KDE));
    if (!kde) {
        return NULL;
    }
    
    kde->n_dims = n_dims;
    kde->bandwidth = KDE_DEFAULT_BANDWIDTH;
    kde->bandwidths = NULL;  /* Will be allocated in kde_fit */
    kde->data = NULL;
    kde->n_points = 0;
    
    return kde;
}

/**
 * @brief Destroy a KDE estimator
 * 
 * @param kde KDE to destroy
 */
void kde_destroy(KDE *kde) {
    if (kde) {
        free(kde->data);
        free(kde->bandwidths);
        free(kde);
    }
}

/**
 * @brief Compute optimal bandwidth using Silverman's rule of thumb
 * 
 * h = 0.9 * min(std, IQR/1.34) * n^(-1/5)
 * 
 * @param data Data array
 * @param n Number of points
 * @return Optimal bandwidth
 */
double kde_bandwidth_silverman(const double *data, uint32_t n) {
    if (n < 2) {
        return KDE_DEFAULT_BANDWIDTH;
    }
    
    /* Compute mean */
    double sum = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        sum += data[i];
    }
    double mean_val = sum / n;
    
    /* Compute standard deviation */
    double var_sum = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        double diff = data[i] - mean_val;
        var_sum += diff * diff;
    }
    double std = sqrt(var_sum / (n - 1));
    
    /* Copy data for IQR calculation */
    double *sorted = (double *)malloc(n * sizeof(double));
    if (!sorted) {
        return 0.9 * std * pow(n, -0.2);
    }
    memcpy(sorted, data, n * sizeof(double));
    
    /* Sort for IQR */
    for (uint32_t i = 0; i < n - 1; i++) {
        for (uint32_t j = i + 1; j < n; j++) {
            if (sorted[j] < sorted[i]) {
                double tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }
    
    /* Compute IQR */
    uint32_t q1_idx = n / 4;
    uint32_t q3_idx = 3 * n / 4;
    double iqr = sorted[q3_idx] - sorted[q1_idx];
    
    free(sorted);
    
    /* Silverman's rule */
    double scale = fmin(std, iqr / 1.34);
    if (scale < 1e-10) {
        scale = std;
    }
    
    return 0.9 * scale * pow(n, -0.2);
}

/**
 * @brief Fit the KDE to data with per-dimension bandwidths
 * 
 * Uses Silverman's rule for each dimension, but enforces minimum
 * bandwidths for sky position dimensions to smooth over survey
 * artifacts like DESI's hexagonal tile pattern (~1.6 deg tiles).
 * 
 * Dimension layout expected: [phi, lambda, mu_perp, ...]
 * - Dimensions 0,1 (sky position): min 2.0 deg bandwidth
 * - Other dimensions: Silverman's rule only
 * 
 * @param kde KDE to fit
 * @param data 2D array of data points (n_points x n_dims)
 * @param n_points Number of data points
 * @return 0 on success, -1 on error
 */
int kde_fit(KDE *kde, double **data, uint32_t n_points) {
    if (!kde || !data || n_points == 0) {
        return -1;
    }
    
    /* Allocate flattened data array */
    kde->data = (double *)malloc(n_points * kde->n_dims * sizeof(double));
    if (!kde->data) {
        return -1;
    }
    
    /* Allocate per-dimension bandwidths */
    kde->bandwidths = (double *)malloc(kde->n_dims * sizeof(double));
    if (!kde->bandwidths) {
        free(kde->data);
        kde->data = NULL;
        return -1;
    }
    
    /* Copy and flatten data */
    for (uint32_t i = 0; i < n_points; i++) {
        for (uint32_t d = 0; d < kde->n_dims; d++) {
            kde->data[i * kde->n_dims + d] = data[i][d];
        }
    }
    
    kde->n_points = n_points;
    
    /* Compute per-dimension bandwidths with minimum thresholds */
    double log_h_sum = 0.0;
    for (uint32_t d = 0; d < kde->n_dims; d++) {
        double *col = (double *)malloc(n_points * sizeof(double));
        if (!col) {
            kde->bandwidths[d] = KDE_DEFAULT_BANDWIDTH;
            log_h_sum += log(KDE_DEFAULT_BANDWIDTH);
            continue;
        }
        
        for (uint32_t i = 0; i < n_points; i++) {
            col[i] = data[i][d];
        }
        
        double h_d = kde_bandwidth_silverman(col, n_points);
        
        /* Apply minimum bandwidth for sky position dimensions (0 and 1)
         * to smooth over DESI tile pattern (~1.6 deg tiles, ~1.4 deg spacing).
         * Using 2.0 deg ensures we don't pick up tile edges as real structure.
         */
        if (d < 2) {
            /* Sky position dimensions: phi and lambda */
            if (h_d < MIN_BANDWIDTH_SKY_DEG) {
                h_d = MIN_BANDWIDTH_SKY_DEG;
            }
        }
        /* PM and other dimensions: use Silverman's rule as-is */
        
        kde->bandwidths[d] = h_d;
        log_h_sum += log(h_d);
        
        free(col);
    }
    
    /* Also store geometric mean for backward compatibility */
    kde->bandwidth = exp(log_h_sum / kde->n_dims);
    
    return 0;
}

/**
 * @brief Evaluate KDE at a point using Gaussian kernel with per-dimension bandwidths
 * 
 * Uses a product of 1D Gaussian kernels, each with its own bandwidth.
 * This allows different smoothing scales for sky position (large, to smooth
 * over survey artifacts) vs proper motion (data-driven).
 * 
 * @param kde Fitted KDE
 * @param point Point to evaluate at (n_dims values)
 * @return Estimated density
 */
double kde_evaluate(const KDE *kde, const double *point) {
    if (!kde || !kde->data || !point || kde->n_points == 0) {
        return 0.0;
    }
    
    /* Use per-dimension bandwidths if available */
    if (kde->bandwidths) {
        /* Product kernel: K(x) = prod_d K_d(x_d / h_d) / h_d */
        double log_norm = 0.0;
        for (uint32_t d = 0; d < kde->n_dims; d++) {
            double h_d = kde->bandwidths[d];
            log_norm += -0.5 * log(2 * M_PI) - log(h_d);
        }
        double norm = exp(log_norm);
        
        double sum = 0.0;
        for (uint32_t i = 0; i < kde->n_points; i++) {
            double log_kernel = 0.0;
            
            for (uint32_t d = 0; d < kde->n_dims; d++) {
                double h_d = kde->bandwidths[d];
                double diff = point[d] - kde->data[i * kde->n_dims + d];
                double z = diff / h_d;
                log_kernel += -0.5 * z * z;
            }
            
            sum += exp(log_kernel);
        }
        
        return norm * sum / kde->n_points;
    }
    
    /* Fallback: use scalar bandwidth (backward compatibility) */
    double h = kde->bandwidth;
    double h_sq = h * h;
    double norm = pow(2 * M_PI * h_sq, -0.5 * kde->n_dims);
    
    double sum = 0.0;
    
    for (uint32_t i = 0; i < kde->n_points; i++) {
        double dist_sq = 0.0;
        
        for (uint32_t d = 0; d < kde->n_dims; d++) {
            double diff = point[d] - kde->data[i * kde->n_dims + d];
            dist_sq += diff * diff;
        }
        
        /* Gaussian kernel */
        sum += exp(-0.5 * dist_sq / h_sq);
    }
    
    return norm * sum / kde->n_points;
}

/* ============================================================================
 * K-Nearest Neighbors Density Estimation
 * ============================================================================ */

/**
 * @brief Compare function for sorting distances
 */
static int compare_doubles(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/**
 * @brief Compute KNN density estimate at a point using kd-tree
 * 
 * @param point Point to evaluate (3D: x, y, z)
 * @param tree Pre-built kd-tree
 * @param k Number of neighbors
 * @param n_dims Number of dimensions (must be 3 for kd-tree)
 * @return KNN density estimate
 */
static double knn_density_kdtree(const double *point, const KDTree3D *tree,
                                  int k, uint32_t n_dims) {
    if (tree->n_points < k) {
        k = tree->n_points;
    }
    if (k < 1) return 0.0;
    
    /* Query k nearest neighbors */
    int *neighbor_idx = (int *)malloc(k * sizeof(int));
    float *neighbor_dist_sq = (float *)malloc(k * sizeof(float));
    if (!neighbor_idx || !neighbor_dist_sq) {
        free(neighbor_idx);
        free(neighbor_dist_sq);
        return 0.0;
    }
    
    int found = kdtree3d_knearest(tree, 
                                  (float)point[0], (float)point[1], (float)point[2],
                                  k, neighbor_idx, neighbor_dist_sq);
    
    if (found < k) k = found;
    if (k < 1) {
        free(neighbor_idx);
        free(neighbor_dist_sq);
        return 0.0;
    }
    
    /* Get k-th nearest neighbor distance */
    double r_k = sqrt((double)neighbor_dist_sq[k - 1]);
    if (r_k < 1e-10) {
        r_k = 1e-10;  /* Avoid division by zero */
    }
    
    free(neighbor_idx);
    free(neighbor_dist_sq);
    
    /* Volume of n-dimensional hypersphere: V = pi^(d/2) * r^d / Gamma(d/2 + 1) */
    double log_vol = (n_dims / 2.0) * log(M_PI) + n_dims * log(r_k);
    /* Approximate Gamma(d/2 + 1) using Stirling */
    double half_d = n_dims / 2.0;
    double log_gamma = half_d * log(half_d) - half_d + 0.5 * log(2 * M_PI / half_d);
    log_vol -= log_gamma;
    
    double vol = exp(log_vol);
    
    /* Density = k / (n * V) */
    return (double)k / (tree->n_points * vol);
}

/* ============================================================================
 * Anomaly Score Computation
 * ============================================================================ */

/**
 * @brief Extract feature vector from a star
 * 
 * Features: (phi, lambda, mu_perp) - position and orthogonal proper motion only
 * No isochrone/photometry features (g_mag, bp_rp removed)
 * 
 * @param star Star to extract features from
 * @param features Output array (must be pre-allocated)
 * @param use_mu_lambda If true, exclude mu_lambda from features (it defines SR)
 */
static void extract_features(const Star *star, double *features, bool use_mu_lambda) {
    features[0] = star->phi;
    features[1] = star->lambda;
    
    if (use_mu_lambda) {
        /* mu_lambda defines SR, so use mu_phi as anomaly feature */
        features[2] = star->mu_phi;
    } else {
        /* mu_phi defines SR, so use mu_lambda as anomaly feature */
        features[2] = star->mu_lambda;
    }
    /* No photometry features - algorithm uses only positions and proper motions */
}

#define N_FEATURES 3

/**
 * @brief Compute anomaly scores using KDE
 * 
 * Trains KDE on sideband and evaluates density ratio for signal region stars.
 * Uses per-dimension bandwidths with minimum scale for sky positions to
 * smooth over survey artifacts.
 * 
 * @param patch Patch containing stars
 * @param sr Signal region definition
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int compute_anomaly_scores_kde(Patch *patch, SignalRegion *sr, const Config *cfg) {
    if (!patch || !sr || sr->n_sb_stars < 100) {
        return -1;
    }
    
    /* Get minimum sky bandwidth from config (default 2.0 deg) */
    double min_sky_bw = (cfg && cfg->min_sky_bandwidth > 0) ? 
                         cfg->min_sky_bandwidth : MIN_BANDWIDTH_SKY_DEG;
    
    printf("Computing KDE anomaly scores (min_sky_bandwidth=%.1f deg)...\n", min_sky_bw);
    printf("  SR stars: %u, SB stars: %u\n", sr->n_sr_stars, sr->n_sb_stars);
    
    /* Extract features from sideband */
    double **sb_features = (double **)malloc(sr->n_sb_stars * sizeof(double *));
    if (!sb_features) return -1;
    
    for (uint32_t i = 0; i < sr->n_sb_stars; i++) {
        sb_features[i] = (double *)malloc(N_FEATURES * sizeof(double));
        if (!sb_features[i]) {
            for (uint32_t j = 0; j < i; j++) free(sb_features[j]);
            free(sb_features);
            return -1;
        }
        extract_features(sr->sb_stars[i], sb_features[i], sr->use_mu_lambda);
    }
    
    /* Fit KDE to sideband (background model) */
    KDE *kde_bg = kde_create(N_FEATURES);
    if (!kde_bg) {
        for (uint32_t i = 0; i < sr->n_sb_stars; i++) free(sb_features[i]);
        free(sb_features);
        return -1;
    }
    
    if (kde_fit(kde_bg, sb_features, sr->n_sb_stars) != 0) {
        kde_destroy(kde_bg);
        for (uint32_t i = 0; i < sr->n_sb_stars; i++) free(sb_features[i]);
        free(sb_features);
        return -1;
    }
    
    /* Extract features from signal region */
    double **sr_features = (double **)malloc(sr->n_sr_stars * sizeof(double *));
    if (!sr_features) {
        kde_destroy(kde_bg);
        for (uint32_t i = 0; i < sr->n_sb_stars; i++) free(sb_features[i]);
        free(sb_features);
        return -1;
    }
    
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        sr_features[i] = (double *)malloc(N_FEATURES * sizeof(double));
        if (!sr_features[i]) {
            for (uint32_t j = 0; j < i; j++) free(sr_features[j]);
            free(sr_features);
            kde_destroy(kde_bg);
            for (uint32_t k = 0; k < sr->n_sb_stars; k++) free(sb_features[k]);
            free(sb_features);
            return -1;
        }
        extract_features(sr->sr_stars[i], sr_features[i], sr->use_mu_lambda);
    }
    
    /* Fit KDE to signal region (data model) */
    KDE *kde_data = kde_create(N_FEATURES);
    if (!kde_data || kde_fit(kde_data, sr_features, sr->n_sr_stars) != 0) {
        if (kde_data) kde_destroy(kde_data);
        kde_destroy(kde_bg);
        for (uint32_t i = 0; i < sr->n_sr_stars; i++) free(sr_features[i]);
        free(sr_features);
        for (uint32_t i = 0; i < sr->n_sb_stars; i++) free(sb_features[i]);
        free(sb_features);
        return -1;
    }
    
    /* Compute anomaly scores for signal region stars */
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        double p_data = kde_evaluate(kde_data, sr_features[i]);
        double p_bg = kde_evaluate(kde_bg, sr_features[i]);
        
        /* Anomaly score R(x) = p_data(x) / p_bg(x) */
        if (p_bg > 1e-300) {
            sr->sr_stars[i]->anomaly_score = p_data / p_bg;
        } else {
            sr->sr_stars[i]->anomaly_score = (p_data > 1e-300) ? 1e10 : 1.0;
        }
    }
    
    /* Cleanup */
    kde_destroy(kde_data);
    kde_destroy(kde_bg);
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) free(sr_features[i]);
    free(sr_features);
    for (uint32_t i = 0; i < sr->n_sb_stars; i++) free(sb_features[i]);
    free(sb_features);
    
    return 0;
}

/**
 * @brief Compute anomaly scores using KNN density estimation with kd-tree
 * 
 * Uses a kd-tree for O(log N) neighbor queries instead of O(N) brute force.
 * Computes local density and flags overdense stars.
 * 
 * IMPORTANT: Uses standardization with minimum scale for spatial dimensions
 * to smooth over survey artifacts. The effective smoothing scale for sky
 * positions is max(data_std, MIN_BANDWIDTH_SKY_DEG), ensuring we don't
 * pick up DESI tile edges as real overdensities.
 * 
 * @param patch Patch containing stars
 * @param sr Signal region definition
 * @param k Number of neighbors for KNN
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int compute_anomaly_scores_knn(Patch *patch, SignalRegion *sr, int k, const Config *cfg) {
    (void)patch;
    
    if (!sr || sr->n_sr_stars < 10) {
        return -1;
    }
    
    /* Adjust k if necessary */
    if (k > (int)sr->n_sr_stars / 2) {
        k = sr->n_sr_stars / 2;
    }
    if (k < 3) k = 3;
    
    /* Determine which dimensions to use based on config */
    bool use_dist = cfg && cfg->use_distance_knn;
    bool use_rv = cfg && cfg->use_rv_knn;
    
    /* Allocate arrays for kd-tree coordinates */
    float *x = (float *)malloc(sr->n_sr_stars * sizeof(float));
    float *y = (float *)malloc(sr->n_sr_stars * sizeof(float));
    float *z = (float *)malloc(sr->n_sr_stars * sizeof(float));
    
    if (!x || !y || !z) {
        free(x); free(y); free(z);
        return -1;
    }
    
    /* 
     * Variance-based standardization for unit-less KNN
     * Each dimension is standardized: (value - mean) / std
     * This ensures all dimensions contribute equally regardless of units
     */
    
    /* First pass: compute means for each dimension */
    double mean_phi = 0.0, mean_lambda = 0.0, mean_dist = 0.0, mean_rv = 0.0;
    int n_valid_dist = 0, n_valid_rv = 0;
    
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        mean_phi += sr->sr_stars[i]->phi;
        mean_lambda += sr->sr_stars[i]->lambda;
        
        if (use_dist && sr->sr_stars[i]->distance > 0.1 && sr->sr_stars[i]->distance < 500.0) {
            mean_dist += sr->sr_stars[i]->distance;
            n_valid_dist++;
        }
        if (use_rv && fabs(sr->sr_stars[i]->rv) > 0.001 && fabs(sr->sr_stars[i]->rv) < 1000.0) {
            mean_rv += sr->sr_stars[i]->rv;
            n_valid_rv++;
        }
    }
    mean_phi /= sr->n_sr_stars;
    mean_lambda /= sr->n_sr_stars;
    if (n_valid_dist > 0) mean_dist /= n_valid_dist;
    if (n_valid_rv > 0) mean_rv /= n_valid_rv;
    
    /* Second pass: compute variances */
    double var_phi = 0.0, var_lambda = 0.0, var_dist = 0.0, var_rv = 0.0;
    
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        var_phi += pow(sr->sr_stars[i]->phi - mean_phi, 2);
        var_lambda += pow(sr->sr_stars[i]->lambda - mean_lambda, 2);
        
        if (use_dist && sr->sr_stars[i]->distance > 0.1 && sr->sr_stars[i]->distance < 500.0) {
            var_dist += pow(sr->sr_stars[i]->distance - mean_dist, 2);
        }
        if (use_rv && fabs(sr->sr_stars[i]->rv) > 0.001 && fabs(sr->sr_stars[i]->rv) < 1000.0) {
            var_rv += pow(sr->sr_stars[i]->rv - mean_rv, 2);
        }
    }
    var_phi /= sr->n_sr_stars;
    var_lambda /= sr->n_sr_stars;
    if (n_valid_dist > 1) var_dist /= n_valid_dist;
    if (n_valid_rv > 1) var_rv /= n_valid_rv;
    
    /* Standard deviations (avoid division by zero)
     * 
     * IMPORTANT: Use minimum smoothing scale for spatial dimensions to smooth
     * over survey artifacts like DESI's hexagonal tile pattern (~1.6 deg tiles).
     * This prevents tile edges from being flagged as real overdensities.
     */
    double std_phi = sqrt(var_phi > 1e-10 ? var_phi : 1e-10);
    double std_lambda = sqrt(var_lambda > 1e-10 ? var_lambda : 1e-10);
    
    /* Get minimum sky bandwidth from config (default 2.0 deg) */
    double min_sky_bw = (cfg && cfg->min_sky_bandwidth > 0) ? 
                         cfg->min_sky_bandwidth : MIN_BANDWIDTH_SKY_DEG;
    
    /* Enforce minimum smoothing scale for sky positions */
    if (std_phi < min_sky_bw) {
        std_phi = min_sky_bw;
    }
    if (std_lambda < min_sky_bw) {
        std_lambda = min_sky_bw;
    }
    
    double std_dist = sqrt(var_dist > 1e-10 ? var_dist : 1.0);  /* Default 1 if no variance */
    double std_rv = sqrt(var_rv > 1e-10 ? var_rv : 1.0);
    
    /* Third pass: fill standardized coordinates
     * x = standardized phi
     * y = standardized lambda  
     * z = standardized distance (if enabled) or standardized RV (if only RV enabled)
     *     or combination if both enabled (we'll combine dist+rv into z using weighted sum)
     */
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        x[i] = (float)((sr->sr_stars[i]->phi - mean_phi) / std_phi);
        y[i] = (float)((sr->sr_stars[i]->lambda - mean_lambda) / std_lambda);
        
        if (use_dist && use_rv) {
            /* Both distance and RV: combine into z as weighted average of standardized values */
            double std_d = 0.0, std_r = 0.0;
            double dist = sr->sr_stars[i]->distance;
            double rv = sr->sr_stars[i]->rv;
            
            /* Standardize distance (or use 0 if invalid - neutral) */
            if (dist > 0.1 && dist < 500.0) {
                std_d = (dist - mean_dist) / std_dist;
            }
            
            /* Standardize RV (or use 0 if invalid - neutral) */
            if (fabs(rv) > 0.001 && fabs(rv) < 1000.0) {
                std_r = (rv - mean_rv) / std_rv;
            }
            
            /* Combine: use Euclidean combination so both contribute */
            z[i] = (float)(sqrt(std_d * std_d + std_r * std_r) / sqrt(2.0));
            
        } else if (use_dist) {
            /* Only distance */
            double dist = sr->sr_stars[i]->distance;
            if (dist > 0.1 && dist < 500.0) {
                z[i] = (float)((dist - mean_dist) / std_dist);
            } else {
                z[i] = 0.0f;  /* Neutral if invalid */
            }
            
        } else if (use_rv) {
            /* Only RV */
            double rv = sr->sr_stars[i]->rv;
            if (fabs(rv) > 0.001 && fabs(rv) < 1000.0) {
                z[i] = (float)((rv - mean_rv) / std_rv);
            } else {
                z[i] = 0.0f;  /* Neutral if invalid */
            }
            
        } else {
            /* Neither: 2D mode (original behavior) */
            z[i] = 0.0f;
        }
    }
    
    /* Build the kd-tree with standardized coordinates */
    KDTree3D tree;
    kdtree3d_init(&tree, 1.0f, 1.0f, 1.0f);  /* All dimensions equally weighted now */
    if (kdtree3d_build(x, y, z, (int)sr->n_sr_stars, &tree) != 0) {
        free(x); free(y); free(z);
        return -1;
    }
    
    /* Allocate result arrays */
    double *r_k = (double *)malloc(sr->n_sr_stars * sizeof(double));
    int *neighbor_idx = (int *)malloc((k + 1) * sizeof(int));
    float *neighbor_dist_sq = (float *)malloc((k + 1) * sizeof(float));
    
    if (!r_k || !neighbor_idx || !neighbor_dist_sq) {
        free(r_k); free(neighbor_idx); free(neighbor_dist_sq);
        free(x); free(y); free(z);
        kdtree3d_free(&tree);
        return -1;
    }
    
    /* Query k+1 nearest neighbors for each star (including self) */
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        int found = kdtree3d_knearest(&tree, x[i], y[i], z[i], 
                                       k + 1, neighbor_idx, neighbor_dist_sq);
        
        /* The first neighbor is the point itself (distance 0), so use k-th after that */
        int kth_idx = (found > k) ? k : found - 1;
        if (kth_idx < 0) kth_idx = 0;
        
        r_k[i] = sqrt((double)neighbor_dist_sq[kth_idx]);
        if (r_k[i] < 0.01) r_k[i] = 0.01;  /* Avoid tiny distances */
    }
    
    /* Compute statistics of r_k */
    double sum_r = 0.0, sum_r2 = 0.0;
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        sum_r += r_k[i];
        sum_r2 += r_k[i] * r_k[i];
    }
    double mean_r = sum_r / sr->n_sr_stars;
    double var_r = sum_r2 / sr->n_sr_stars - mean_r * mean_r;
    double std_r = sqrt(var_r > 0 ? var_r : 0.01);
    
    /* Assign anomaly scores: stars with small r_k (high local density) get high scores */
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        /* Anomaly score = how much denser than average */
        /* Small r_k means high density */
        double z_score = (mean_r - r_k[i]) / std_r;  /* Positive z = overdense */
        sr->sr_stars[i]->anomaly_score = z_score;
    }
    
    /* Cleanup */
    free(r_k);
    free(neighbor_idx);
    free(neighbor_dist_sq);
    free(x); free(y); free(z);
    kdtree3d_free(&tree);
    
    return 0;
}

/**
 * @brief Compute anomaly scores using histogram-based density
 * 
 * Simple and fast but less accurate than KDE/KNN.
 * 
 * @param patch Patch containing stars
 * @param sr Signal region definition
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int compute_anomaly_scores_histogram(Patch *patch, SignalRegion *sr, const Config *cfg) {
    (void)patch;  /* Currently unused */
    (void)cfg;    /* Currently unused */
    
    if (!sr || sr->n_sb_stars < 100) {
        return -1;
    }
    
    printf("Computing histogram-based anomaly scores...\n");
    
    /* Use a 5D histogram with coarse bins */
    #define N_BINS 10
    #define TOTAL_BINS (N_BINS * N_BINS * N_BINS * N_BINS * N_BINS)
    
    /* Find data range from sideband */
    double mins[N_FEATURES], maxs[N_FEATURES];
    for (int d = 0; d < N_FEATURES; d++) {
        mins[d] = DBL_MAX;
        maxs[d] = -DBL_MAX;
    }
    
    for (uint32_t i = 0; i < sr->n_sb_stars; i++) {
        double features[N_FEATURES];
        extract_features(sr->sb_stars[i], features, sr->use_mu_lambda);
        for (int d = 0; d < N_FEATURES; d++) {
            if (features[d] < mins[d]) mins[d] = features[d];
            if (features[d] > maxs[d]) maxs[d] = features[d];
        }
    }
    
    /* Add small margin */
    for (int d = 0; d < N_FEATURES; d++) {
        double range = maxs[d] - mins[d];
        if (range < 1e-10) range = 1.0;
        mins[d] -= 0.1 * range;
        maxs[d] += 0.1 * range;
    }
    
    /* Compute bin widths */
    double bin_widths[N_FEATURES];
    for (int d = 0; d < N_FEATURES; d++) {
        bin_widths[d] = (maxs[d] - mins[d]) / N_BINS;
    }
    
    /* Build histogram for sideband */
    int *hist_sb = (int *)calloc(TOTAL_BINS, sizeof(int));
    if (!hist_sb) return -1;
    
    for (uint32_t i = 0; i < sr->n_sb_stars; i++) {
        double features[N_FEATURES];
        extract_features(sr->sb_stars[i], features, sr->use_mu_lambda);
        
        int idx = 0;
        int multiplier = 1;
        for (int d = 0; d < N_FEATURES; d++) {
            int bin = (int)((features[d] - mins[d]) / bin_widths[d]);
            if (bin < 0) bin = 0;
            if (bin >= N_BINS) bin = N_BINS - 1;
            idx += bin * multiplier;
            multiplier *= N_BINS;
        }
        
        hist_sb[idx]++;
    }
    
    /* Build histogram for signal region */
    int *hist_sr = (int *)calloc(TOTAL_BINS, sizeof(int));
    if (!hist_sr) {
        free(hist_sb);
        return -1;
    }
    
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        double features[N_FEATURES];
        extract_features(sr->sr_stars[i], features, sr->use_mu_lambda);
        
        int idx = 0;
        int multiplier = 1;
        for (int d = 0; d < N_FEATURES; d++) {
            int bin = (int)((features[d] - mins[d]) / bin_widths[d]);
            if (bin < 0) bin = 0;
            if (bin >= N_BINS) bin = N_BINS - 1;
            idx += bin * multiplier;
            multiplier *= N_BINS;
        }
        
        hist_sr[idx]++;
    }
    
    /* Compute anomaly scores */
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        double features[N_FEATURES];
        extract_features(sr->sr_stars[i], features, sr->use_mu_lambda);
        
        int idx = 0;
        int multiplier = 1;
        for (int d = 0; d < N_FEATURES; d++) {
            int bin = (int)((features[d] - mins[d]) / bin_widths[d]);
            if (bin < 0) bin = 0;
            if (bin >= N_BINS) bin = N_BINS - 1;
            idx += bin * multiplier;
            multiplier *= N_BINS;
        }
        
        double p_data = (double)(hist_sr[idx] + 1) / (sr->n_sr_stars + TOTAL_BINS);
        double p_bg = (double)(hist_sb[idx] + 1) / (sr->n_sb_stars + TOTAL_BINS);
        
        sr->sr_stars[i]->anomaly_score = p_data / p_bg;
    }
    
    free(hist_sr);
    free(hist_sb);
    
    #undef N_BINS
    #undef TOTAL_BINS
    
    return 0;
}

/* ============================================================================
 * Signal Region and ROI Management
 * ============================================================================ */

/**
 * @brief Create a signal region
 * 
 * IMPORTANT: Uses equatorial proper motions (pmra, pmdec) instead of
 * patch-local (mu_phi, mu_lambda) because stream stars have coherent
 * proper motions in equatorial coordinates.
 * 
 * @param patch Parent patch
 * @param mu_min Minimum proper motion for SR
 * @param mu_max Maximum proper motion for SR
 * @param use_mu_lambda If true, use pmdec; otherwise pmra (equatorial frame)
 * @return Newly allocated SignalRegion, or NULL on error
 */
SignalRegion *sr_create(Patch *patch, double mu_min, double mu_max, bool use_mu_lambda) {
    if (!patch) return NULL;
    
    SignalRegion *sr = (SignalRegion *)calloc(1, sizeof(SignalRegion));
    if (!sr) return NULL;
    
    sr->mu_min = mu_min;
    sr->mu_max = mu_max;
    sr->use_mu_lambda = use_mu_lambda;
    
    /* Count stars in SR and SB */
    /* Use equatorial proper motions (pmra, pmdec) for coherent stream detection */
    uint32_t n_sr = 0, n_sb = 0;
    for (uint32_t i = 0; i < patch->n_stars; i++) {
        /* use_mu_lambda: true=pmdec, false=pmra (in equatorial frame) */
        double mu = use_mu_lambda ? patch->stars[i].pmdec : patch->stars[i].pmra;
        if (mu >= mu_min && mu < mu_max) {
            n_sr++;
        } else {
            n_sb++;
        }
    }
    
    /* Allocate arrays */
    sr->sr_stars = (Star **)malloc(n_sr * sizeof(Star *));
    sr->sb_stars = (Star **)malloc(n_sb * sizeof(Star *));
    if (!sr->sr_stars || !sr->sb_stars) {
        free(sr->sr_stars);
        free(sr->sb_stars);
        free(sr);
        return NULL;
    }
    
    /* Assign stars */
    sr->n_sr_stars = 0;
    sr->n_sb_stars = 0;
    for (uint32_t i = 0; i < patch->n_stars; i++) {
        double mu = use_mu_lambda ? patch->stars[i].pmdec : patch->stars[i].pmra;
        if (mu >= mu_min && mu < mu_max) {
            sr->sr_stars[sr->n_sr_stars++] = &patch->stars[i];
        } else {
            sr->sb_stars[sr->n_sb_stars++] = &patch->stars[i];
        }
    }
    
    return sr;
}

/**
 * @brief Destroy a signal region
 * 
 * @param sr SignalRegion to destroy
 */
void sr_destroy(SignalRegion *sr) {
    if (sr) {
        free(sr->sr_stars);
        free(sr->sb_stars);
        free(sr);
    }
}

/**
 * @brief Create a region of interest within a signal region
 * 
 * IMPORTANT: Uses equatorial proper motions (pmra, pmdec) for coherent
 * stream detection. If SR uses pmdec, ROI uses pmra and vice versa.
 * 
 * @param sr Parent signal region
 * @param mu_perp_min Minimum orthogonal proper motion
 * @param mu_perp_max Maximum orthogonal proper motion
 * @return Newly allocated ROI, or NULL on error
 */
ROI *roi_create(SignalRegion *sr, double mu_perp_min, double mu_perp_max) {
    if (!sr) return NULL;
    
    ROI *roi = (ROI *)calloc(1, sizeof(ROI));
    if (!roi) return NULL;
    
    roi->mu_perp_min = mu_perp_min;
    roi->mu_perp_max = mu_perp_max;
    
    /* Count stars in ROI */
    /* Use equatorial proper motions: if SR uses pmdec, ROI uses pmra and vice versa */
    uint32_t n_roi = 0;
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        Star *star = sr->sr_stars[i];
        if (!star->passes_fiducial) continue;
        
        /* use_mu_lambda: if SR used pmdec, ROI uses pmra (and vice versa) */
        double mu_perp = sr->use_mu_lambda ? star->pmra : star->pmdec;
        if (mu_perp >= mu_perp_min && mu_perp < mu_perp_max) {
            n_roi++;
        }
    }
    
    if (n_roi < 5) {
        /* Not enough stars in ROI */
        free(roi);
        return NULL;
    }
    
    /* Allocate arrays */
    roi->roi_stars = (Star **)malloc(n_roi * sizeof(Star *));
    if (!roi->roi_stars) {
        free(roi);
        return NULL;
    }
    
    /* Assign stars */
    /* Use equatorial proper motions for orthogonal direction */
    roi->n_roi_stars = 0;
    for (uint32_t i = 0; i < sr->n_sr_stars; i++) {
        Star *star = sr->sr_stars[i];
        if (!star->passes_fiducial) continue;
        
        double mu_perp = sr->use_mu_lambda ? star->pmra : star->pmdec;
        if (mu_perp >= mu_perp_min && mu_perp < mu_perp_max) {
            roi->roi_stars[roi->n_roi_stars++] = star;
        }
    }
    
    /* Select ALL stars in the ROI for line detection
     * For stream detection, we need to find linear structures, not clusters.
     * Using all stars in the proper motion window gives the best chance
     * of detecting streams as lines in position space.
     * 
     * If there are too many stars (>500), we subsample to keep computation fast.
     */
    int n_select = roi->n_roi_stars;
    int max_stars_for_hough = 500;  /* Max stars for Hough transform */
    
    /* Minimum 5 stars needed for meaningful line detection */
    if (n_select < 5) {
        free(roi->roi_stars);
        free(roi);
        return NULL;
    }
    
    roi->selected_stars = (Star **)malloc(n_select * sizeof(Star *));
    if (!roi->selected_stars) {
        free(roi->roi_stars);
        free(roi);
        return NULL;
    }
    
    if (n_select <= max_stars_for_hough) {
        /* Use all stars in the ROI */
        for (uint32_t i = 0; i < roi->n_roi_stars; i++) {
            roi->selected_stars[i] = roi->roi_stars[i];
        }
    } else {
        /* Subsample: take every n-th star to get ~max_stars_for_hough */
        int step = n_select / max_stars_for_hough;
        if (step < 1) step = 1;
        n_select = 0;
        for (uint32_t i = 0; i < roi->n_roi_stars && n_select < max_stars_for_hough; i += step) {
            roi->selected_stars[n_select++] = roi->roi_stars[i];
        }
    }
    
    roi->n_selected = n_select;
    
    return roi;
}

/**
 * @brief Destroy a region of interest
 * 
 * @param roi ROI to destroy
 */
void roi_destroy(ROI *roi) {
    if (roi) {
        free(roi->roi_stars);
        free(roi->selected_stars);
        free(roi);
    }
}

