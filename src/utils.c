/**
 * @file utils.c
 * @brief Utility functions for stream detection
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stream_detect.h"

/* ============================================================================
 * Memory Management
 * ============================================================================ */

/**
 * @brief Safe malloc with error checking
 */
void *safe_malloc(size_t size) {
    void *ptr = malloc(size);
    if (!ptr && size > 0) {
        fprintf(stderr, "Error: Memory allocation failed (size=%zu)\n", size);
        exit(1);
    }
    return ptr;
}

/**
 * @brief Safe calloc with error checking
 */
void *safe_calloc(size_t nmemb, size_t size) {
    void *ptr = calloc(nmemb, size);
    if (!ptr && nmemb > 0 && size > 0) {
        fprintf(stderr, "Error: Memory allocation failed (nmemb=%zu, size=%zu)\n", nmemb, size);
        exit(1);
    }
    return ptr;
}

/**
 * @brief Safe realloc with error checking
 */
void *safe_realloc(void *ptr, size_t size) {
    void *new_ptr = realloc(ptr, size);
    if (!new_ptr && size > 0) {
        fprintf(stderr, "Error: Memory reallocation failed (size=%zu)\n", size);
        exit(1);
    }
    return new_ptr;
}

/* ============================================================================
 * Statistical Functions
 * ============================================================================ */

/**
 * @brief Compute mean of array
 */
double mean(const double *data, uint32_t n) {
    if (n == 0) return 0.0;
    
    double sum = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / n;
}

/**
 * @brief Compute variance of array
 */
double variance(const double *data, uint32_t n) {
    if (n < 2) return 0.0;
    
    double m = mean(data, n);
    double sum_sq = 0.0;
    
    for (uint32_t i = 0; i < n; i++) {
        double diff = data[i] - m;
        sum_sq += diff * diff;
    }
    
    return sum_sq / (n - 1);  /* Sample variance */
}

/**
 * @brief Compute standard deviation of array
 */
double std_dev(const double *data, uint32_t n) {
    return sqrt(variance(data, n));
}

/**
 * @brief Compute median of array (modifies array order)
 */
double median(double *data, uint32_t n) {
    if (n == 0) return 0.0;
    if (n == 1) return data[0];
    
    /* Sort using insertion sort (efficient for small arrays) */
    for (uint32_t i = 1; i < n; i++) {
        double key = data[i];
        int j = i - 1;
        while (j >= 0 && data[j] > key) {
            data[j + 1] = data[j];
            j--;
        }
        data[j + 1] = key;
    }
    
    if (n % 2 == 0) {
        return (data[n/2 - 1] + data[n/2]) / 2.0;
    } else {
        return data[n/2];
    }
}

/**
 * @brief Compute percentile of array (modifies array order)
 */
double percentile(double *data, uint32_t n, double p) {
    if (n == 0) return 0.0;
    if (n == 1) return data[0];
    if (p <= 0) return data[0];
    if (p >= 100) return data[n-1];
    
    /* Sort */
    for (uint32_t i = 1; i < n; i++) {
        double key = data[i];
        int j = i - 1;
        while (j >= 0 && data[j] > key) {
            data[j + 1] = data[j];
            j--;
        }
        data[j + 1] = key;
    }
    
    /* Linear interpolation */
    double rank = (p / 100.0) * (n - 1);
    uint32_t lower = (uint32_t)rank;
    uint32_t upper = lower + 1;
    double frac = rank - lower;
    
    if (upper >= n) {
        return data[n - 1];
    }
    
    return data[lower] * (1 - frac) + data[upper] * frac;
}

/* ============================================================================
 * Sorting Comparators
 * ============================================================================ */

/**
 * @brief Compare stars by anomaly score (descending)
 */
int compare_stars_by_anomaly_score(const void *a, const void *b) {
    const Star *sa = *(const Star **)a;
    const Star *sb = *(const Star **)b;
    
    if (sa->anomaly_score > sb->anomaly_score) return -1;
    if (sa->anomaly_score < sb->anomaly_score) return 1;
    return 0;
}

/**
 * @brief Compare stars by magnitude (ascending)
 */
int compare_stars_by_magnitude(const void *a, const void *b) {
    const Star *sa = *(const Star **)a;
    const Star *sb = *(const Star **)b;
    
    if (sa->g_mag < sb->g_mag) return -1;
    if (sa->g_mag > sb->g_mag) return 1;
    return 0;
}

/* ============================================================================
 * Additional Utility Functions
 * ============================================================================ */

/**
 * @brief Normalize angle to [0, 360)
 */
double normalize_angle_360(double angle) {
    while (angle < 0) angle += 360.0;
    while (angle >= 360) angle -= 360.0;
    return angle;
}

/**
 * @brief Normalize angle to [-180, 180)
 */
double normalize_angle_180(double angle) {
    while (angle < -180) angle += 360.0;
    while (angle >= 180) angle -= 360.0;
    return angle;
}

/**
 * @brief Clip value to range [min, max]
 */
double clip(double value, double min_val, double max_val) {
    if (value < min_val) return min_val;
    if (value > max_val) return max_val;
    return value;
}

/**
 * @brief Linear interpolation
 */
double lerp(double a, double b, double t) {
    return a + t * (b - a);
}

/**
 * @brief Check if value is NaN
 */
int is_nan(double x) {
    return x != x;  /* NaN is the only value not equal to itself */
}

/**
 * @brief Check if value is finite (not NaN or Inf)
 */
int is_finite(double x) {
    return !is_nan(x) && x != INFINITY && x != -INFINITY;
}

