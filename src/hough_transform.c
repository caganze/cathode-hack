/**
 * @file hough_transform.c
 * @brief Hough transform for line detection in stellar streams
 * 
 * Implements the Hough transform to find line-like structures
 * among anomalous stars in position space (phi, lambda).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stream_detect.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Hough Accumulator Management
 * ============================================================================ */

/**
 * @brief Create a Hough accumulator
 * 
 * @param n_theta Number of theta bins
 * @param n_rho Number of rho bins
 * @param rho_max Maximum rho value (usually patch radius)
 * @return Newly allocated HoughAccumulator, or NULL on error
 */
HoughAccumulator *hough_create(int n_theta, int n_rho, double rho_max) {
    HoughAccumulator *ha = (HoughAccumulator *)calloc(1, sizeof(HoughAccumulator));
    if (!ha) {
        return NULL;
    }
    
    ha->n_theta = n_theta;
    ha->n_rho = n_rho;
    ha->theta_min = -M_PI;
    ha->theta_max = M_PI;
    ha->rho_max = rho_max;
    
    ha->d_theta = (ha->theta_max - ha->theta_min) / n_theta;
    ha->d_rho = (2.0 * rho_max) / n_rho;  /* rho ranges from -rho_max to +rho_max */
    
    ha->accumulator = (int *)calloc(n_theta * n_rho, sizeof(int));
    if (!ha->accumulator) {
        free(ha);
        return NULL;
    }
    
    return ha;
}

/**
 * @brief Destroy a Hough accumulator
 * 
 * @param ha HoughAccumulator to destroy
 */
void hough_destroy(HoughAccumulator *ha) {
    if (ha) {
        free(ha->accumulator);
        free(ha);
    }
}

/**
 * @brief Clear the Hough accumulator
 * 
 * @param ha HoughAccumulator to clear
 */
void hough_clear(HoughAccumulator *ha) {
    if (ha && ha->accumulator) {
        memset(ha->accumulator, 0, ha->n_theta * ha->n_rho * sizeof(int));
    }
}

/**
 * @brief Add a point to the Hough accumulator
 * 
 * For each possible line angle theta, compute the perpendicular distance rho
 * and increment the corresponding accumulator bin.
 * 
 * @param ha HoughAccumulator
 * @param x X coordinate (phi in degrees)
 * @param y Y coordinate (lambda in degrees)
 */
void hough_add_point(HoughAccumulator *ha, double x, double y) {
    if (!ha) return;
    
    for (int t = 0; t < ha->n_theta; t++) {
        double theta = ha->theta_min + (t + 0.5) * ha->d_theta;
        
        /* rho = x * cos(theta) + y * sin(theta) */
        double rho = x * cos(theta) + y * sin(theta);
        
        /* Convert rho to bin index */
        int r = (int)((rho + ha->rho_max) / ha->d_rho);
        if (r >= 0 && r < ha->n_rho) {
            ha->accumulator[t * ha->n_rho + r]++;
        }
    }
}

/**
 * @brief Find the peak in the Hough accumulator
 * 
 * @param ha HoughAccumulator
 * @param theta Output best-fit theta (radians)
 * @param rho Output best-fit rho (degrees)
 * @param votes Output number of votes for the peak
 * @return 0 on success, -1 on error
 */
int hough_find_peak(const HoughAccumulator *ha, double *theta, double *rho, int *votes) {
    if (!ha || !theta || !rho || !votes) {
        return -1;
    }
    
    int max_votes = 0;
    int max_t = 0, max_r = 0;
    
    for (int t = 0; t < ha->n_theta; t++) {
        for (int r = 0; r < ha->n_rho; r++) {
            int v = ha->accumulator[t * ha->n_rho + r];
            if (v > max_votes) {
                max_votes = v;
                max_t = t;
                max_r = r;
            }
        }
    }
    
    *theta = ha->theta_min + (max_t + 0.5) * ha->d_theta;
    *rho = -ha->rho_max + (max_r + 0.5) * ha->d_rho;
    *votes = max_votes;
    
    return 0;
}

/**
 * @brief Compute the statistical significance of a Hough line
 * 
 * Uses the binomial distribution to estimate the probability of
 * observing at least 'votes' points on a line by chance.
 * 
 * For a line in 2D, we expect approximately n/sqrt(n_bins) points per bin
 * under random uniform distribution.
 * 
 * @param ha HoughAccumulator (for geometry info)
 * @param votes Number of votes for the line
 * @param n_points Total number of points
 * @return Significance in sigma
 */
double hough_line_significance(const HoughAccumulator *ha, int votes, int n_points) {
    if (!ha || n_points <= 0 || votes <= 0) {
        return 0.0;
    }
    
    /* Expected number of votes per rho bin:
     * Each point votes in all theta bins, creating a sinusoid in (theta, rho) space.
     * The expected density per rho bin is approximately n_points / n_rho
     * But the peak accumulation comes from aligned points.
     */
    double expected = (double)n_points / ha->n_rho;
    if (expected < 1.0) expected = 1.0;
    
    /* Standard deviation under Poisson */
    double sigma = sqrt(expected);
    if (sigma < 1.0) sigma = 1.0;
    
    /* Significance: how many sigma above expected */
    double significance = (double)(votes - expected) / sigma;
    
    /* Also consider: what fraction of points are on the line?
     * If more than 30% of points fall on a single line, that's significant */
    double fraction = (double)votes / n_points;
    if (fraction > 0.3) {
        significance = fmax(significance, 5.0);  /* Boost significance */
    }
    if (fraction > 0.5) {
        significance = fmax(significance, 8.0);
    }
    
    return significance;
}

/**
 * @brief Find the best-fit line through stars in an ROI using Hough transform
 * 
 * @param roi ROI containing selected anomalous stars
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int find_line_in_roi(ROI *roi, const Config *cfg) {
    if (!roi || !cfg || roi->n_selected < 3) {
        return -1;
    }
    
    /* Create Hough accumulator with appropriate resolution */
    double rho_max = cfg->patch_radius * 1.5;  /* Allow some margin */
    HoughAccumulator *ha = hough_create(90, 100, rho_max);  /* Coarser bins for small samples */
    if (!ha) {
        return -1;
    }
    
    /* Add all selected stars to accumulator */
    for (uint32_t i = 0; i < roi->n_selected; i++) {
        Star *star = roi->selected_stars[i];
        hough_add_point(ha, star->phi, star->lambda);
    }
    
    /* Find peak */
    double theta, rho;
    int votes;
    if (hough_find_peak(ha, &theta, &rho, &votes) != 0) {
        hough_destroy(ha);
        return -1;
    }
    
    /* Store results */
    roi->theta_line = theta;
    roi->rho_line = rho;
    roi->sigma_line = hough_line_significance(ha, votes, roi->n_selected);
    
    hough_destroy(ha);
    
    return 0;
}

/* ============================================================================
 * Line Star Selection
 * ============================================================================ */

/**
 * @brief Distance from a point to a line defined by Hough parameters
 * 
 * @param x X coordinate
 * @param y Y coordinate
 * @param theta Line angle
 * @param rho Line distance from origin
 * @return Perpendicular distance to line
 */
static double point_to_line_distance(double x, double y, double theta, double rho) {
    return fabs(x * cos(theta) + y * sin(theta) - rho);
}

/**
 * @brief Select stars that lie on the detected line
 * 
 * @param roi ROI with detected line
 * @param line_stars Output array of stars on the line
 * @param n_line_stars Output number of line stars
 * @param max_distance Maximum distance from line (degrees)
 * @return 0 on success, error code otherwise
 */
int select_line_stars(ROI *roi, Star ***line_stars, uint32_t *n_line_stars,
                      double max_distance) {
    if (!roi || !line_stars || !n_line_stars) {
        return -1;
    }
    
    /* Count stars within distance threshold */
    uint32_t count = 0;
    for (uint32_t i = 0; i < roi->n_selected; i++) {
        Star *star = roi->selected_stars[i];
        double dist = point_to_line_distance(star->phi, star->lambda,
                                             roi->theta_line, roi->rho_line);
        if (dist <= max_distance) {
            count++;
        }
    }
    
    if (count == 0) {
        *line_stars = NULL;
        *n_line_stars = 0;
        return 0;
    }
    
    /* Allocate and fill */
    *line_stars = (Star **)malloc(count * sizeof(Star *));
    if (!*line_stars) {
        return -1;
    }
    
    *n_line_stars = 0;
    for (uint32_t i = 0; i < roi->n_selected; i++) {
        Star *star = roi->selected_stars[i];
        double dist = point_to_line_distance(star->phi, star->lambda,
                                             roi->theta_line, roi->rho_line);
        if (dist <= max_distance) {
            (*line_stars)[(*n_line_stars)++] = star;
        }
    }
    
    return 0;
}

/* ============================================================================
 * Line Merging Utilities
 * ============================================================================ */

/**
 * @brief Check if two lines are compatible for merging
 * 
 * Lines are compatible if their angles and rho values are close.
 * 
 * @param theta1 First line angle
 * @param rho1 First line rho
 * @param theta2 Second line angle
 * @param rho2 Second line rho
 * @param theta_thresh Maximum angle difference (radians)
 * @param rho_thresh Maximum rho difference (degrees)
 * @return true if lines are compatible
 */
bool lines_compatible(double theta1, double rho1, double theta2, double rho2,
                     double theta_thresh, double rho_thresh) {
    /* Normalize angles to [-pi, pi] */
    while (theta1 > M_PI) theta1 -= 2 * M_PI;
    while (theta1 < -M_PI) theta1 += 2 * M_PI;
    while (theta2 > M_PI) theta2 -= 2 * M_PI;
    while (theta2 < -M_PI) theta2 += 2 * M_PI;
    
    /* Handle angle wraparound (lines at theta and theta+pi are the same) */
    double d_theta = fabs(theta1 - theta2);
    if (d_theta > M_PI) {
        d_theta = 2 * M_PI - d_theta;
    }
    
    /* If angles differ by ~pi, lines are the same but with opposite rho sign */
    if (fabs(d_theta - M_PI) < theta_thresh) {
        d_theta = fabs(d_theta - M_PI);
        rho2 = -rho2;
    }
    
    double d_rho = fabs(rho1 - rho2);
    
    return (d_theta < theta_thresh) && (d_rho < rho_thresh);
}

/**
 * @brief Compute overlap fraction between two sets of stars
 * 
 * @param stars1 First set of stars
 * @param n1 Number of stars in first set
 * @param stars2 Second set of stars
 * @param n2 Number of stars in second set
 * @return Fraction of stars1 that are in stars2
 */
double star_set_overlap(Star **stars1, uint32_t n1, Star **stars2, uint32_t n2) {
    /* Null and empty checks */
    if (!stars1 || !stars2 || n1 == 0 || n2 == 0) return 0.0;
    
    int overlap = 0;
    for (uint32_t i = 0; i < n1; i++) {
        if (!stars1[i]) continue;  /* Skip null pointers */
        
        for (uint32_t j = 0; j < n2; j++) {
            if (!stars2[j]) continue;  /* Skip null pointers */
            
            /* Compare by position (should be unique enough) */
            if (fabs(stars1[i]->ra - stars2[j]->ra) < 1e-6 &&
                fabs(stars1[i]->dec - stars2[j]->dec) < 1e-6) {
                overlap++;
                break;
            }
        }
    }
    
    return (double)overlap / n1;
}

/* ============================================================================
 * Weighted Line Fitting (for merged protoclusters)
 * ============================================================================ */

/**
 * @brief Compute weighted average line parameters from multiple ROIs
 * 
 * Weights are based on line significance.
 * 
 * @param rois Array of ROIs
 * @param n_rois Number of ROIs
 * @param theta Output weighted average theta
 * @param rho Output weighted average rho
 * @param combined_sig Output combined significance
 * @return 0 on success, error code otherwise
 */
int combine_line_parameters(ROI **rois, int n_rois,
                           double *theta, double *rho, double *combined_sig) {
    if (!rois || !theta || !rho || !combined_sig || n_rois == 0) {
        return -1;
    }
    
    double sum_weight = 0.0;
    double sum_theta = 0.0;
    double sum_rho = 0.0;
    double sum_sig_sq = 0.0;
    
    /* Use significance as weight */
    for (int i = 0; i < n_rois; i++) {
        double w = rois[i]->sigma_line;
        if (w < 0) w = 0;  /* Ignore negative significances */
        
        sum_weight += w;
        sum_theta += w * rois[i]->theta_line;
        sum_rho += w * rois[i]->rho_line;
        sum_sig_sq += rois[i]->sigma_line * rois[i]->sigma_line;
    }
    
    if (sum_weight < 1e-10) {
        /* No valid weights */
        *theta = rois[0]->theta_line;
        *rho = rois[0]->rho_line;
        *combined_sig = rois[0]->sigma_line;
        return 0;
    }
    
    *theta = sum_theta / sum_weight;
    *rho = sum_rho / sum_weight;
    *combined_sig = sqrt(sum_sig_sq);  /* Quadrature sum */
    
    return 0;
}

/* ============================================================================
 * Line Refinement using RANSAC
 * ============================================================================ */

/**
 * @brief Fit a line to points using RANSAC for outlier rejection
 * 
 * @param stars Array of stars
 * @param n_stars Number of stars
 * @param theta Output best-fit theta
 * @param rho Output best-fit rho
 * @param inlier_threshold Maximum distance for inliers
 * @param n_iterations Number of RANSAC iterations
 * @return Number of inliers
 */
int ransac_line_fit(Star **stars, uint32_t n_stars,
                   double *theta, double *rho,
                   double inlier_threshold, int n_iterations) {
    if (!stars || !theta || !rho || n_stars < 2) {
        return 0;
    }
    
    int best_inliers = 0;
    double best_theta = 0.0, best_rho = 0.0;
    
    for (int iter = 0; iter < n_iterations; iter++) {
        /* Randomly select two points */
        int i1 = rand() % n_stars;
        int i2 = rand() % n_stars;
        while (i2 == i1) {
            i2 = rand() % n_stars;
        }
        
        double x1 = stars[i1]->phi;
        double y1 = stars[i1]->lambda;
        double x2 = stars[i2]->phi;
        double y2 = stars[i2]->lambda;
        
        /* Compute line through these two points */
        double dx = x2 - x1;
        double dy = y2 - y1;
        double len = sqrt(dx*dx + dy*dy);
        if (len < 1e-10) continue;
        
        /* Normal to line: (dy, -dx) / len */
        double nx = dy / len;
        double ny = -dx / len;
        
        /* Theta is angle of normal from x-axis */
        double t = atan2(ny, nx);
        /* Rho is distance from origin along normal */
        double r = x1 * nx + y1 * ny;
        
        /* Count inliers */
        int inliers = 0;
        for (uint32_t j = 0; j < n_stars; j++) {
            double dist = point_to_line_distance(stars[j]->phi, stars[j]->lambda, t, r);
            if (dist <= inlier_threshold) {
                inliers++;
            }
        }
        
        if (inliers > best_inliers) {
            best_inliers = inliers;
            best_theta = t;
            best_rho = r;
        }
    }
    
    *theta = best_theta;
    *rho = best_rho;
    
    return best_inliers;
}

