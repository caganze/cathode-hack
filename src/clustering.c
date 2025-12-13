/**
 * @file clustering.c
 * @brief Clustering algorithms for stream detection
 * 
 * Implements the hierarchical clustering strategy from Via Machinae:
 * 1. ROIs -> Protoclusters (within a patch)
 * 2. Protoclusters -> Protostreams (deduplicated within a patch)
 * 3. Protostreams -> Stream Candidates (merged across patches)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stream_detect.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Forward declarations */
extern double star_set_overlap(Star **stars1, uint32_t n1, Star **stars2, uint32_t n2);
extern bool lines_compatible(double theta1, double rho1, double theta2, double rho2,
                            double theta_thresh, double rho_thresh);

/* ============================================================================
 * Union-Find Helper Functions (non-nested for stability)
 * ============================================================================ */

/**
 * @brief Find root with path compression (iterative to avoid stack overflow)
 */
static int uf_find(int *parent, int x) {
    int root = x;
    while (parent[root] != root) {
        root = parent[root];
    }
    /* Path compression */
    while (parent[x] != root) {
        int next = parent[x];
        parent[x] = root;
        x = next;
    }
    return root;
}

/**
 * @brief Union two sets
 */
static void uf_unite(int *parent, int x, int y) {
    int px = uf_find(parent, x);
    int py = uf_find(parent, y);
    if (px != py) {
        parent[px] = py;
    }
}

/* ============================================================================
 * Protocluster Management
 * ============================================================================ */

/**
 * @brief Create a new protocluster
 * 
 * @return Newly allocated Protocluster, or NULL on error
 */
Protocluster *protocluster_create(void) {
    Protocluster *pc = (Protocluster *)calloc(1, sizeof(Protocluster));
    if (!pc) return NULL;
    
    pc->rois = (ROI **)malloc(MAX_ROIS_PER_PATCH * sizeof(ROI *));
    /* Increased allocation: each ROI can have up to 500 stars after changes */
    pc->line_stars = (Star **)malloc(50000 * sizeof(Star *));
    
    if (!pc->rois || !pc->line_stars) {
        free(pc->rois);
        free(pc->line_stars);
        free(pc);
        return NULL;
    }
    
    pc->n_rois = 0;
    pc->n_line_stars = 0;
    
    return pc;
}

/**
 * @brief Destroy a protocluster
 * 
 * @param pc Protocluster to destroy
 */
void protocluster_destroy(Protocluster *pc) {
    if (pc) {
        free(pc->rois);
        free(pc->line_stars);
        free(pc);
    }
}

/**
 * @brief Add an ROI to a protocluster
 * 
 * @param pc Protocluster
 * @param roi ROI to add
 * @return 0 on success, -1 on error
 */
int protocluster_add_roi(Protocluster *pc, ROI *roi) {
    if (!pc || !roi) return -1;
    
    /* Bounds check for ROIs array */
    if (pc->n_rois >= MAX_ROIS_PER_PATCH) return -1;
    
    pc->rois[pc->n_rois++] = roi;
    
    /* Add line stars from this ROI with bounds checking */
    #define MAX_PROTOCLUSTER_STARS 50000
    for (uint32_t i = 0; i < roi->n_selected; i++) {
        /* Bounds check */
        if (pc->n_line_stars >= MAX_PROTOCLUSTER_STARS) break;
        
        /* Check if star already in list */
        bool found = false;
        for (uint32_t j = 0; j < pc->n_line_stars; j++) {
            if (pc->line_stars[j] == roi->selected_stars[i]) {
                found = true;
                break;
            }
        }
        if (!found) {
            pc->line_stars[pc->n_line_stars++] = roi->selected_stars[i];
        }
    }
    #undef MAX_PROTOCLUSTER_STARS
    
    return 0;
}

/**
 * @brief Compute the combined significance of a protocluster
 * 
 * Uses quadrature sum of individual ROI significances.
 * 
 * @param pc Protocluster
 * @return Combined significance
 */
double protocluster_compute_significance(Protocluster *pc) {
    if (!pc || pc->n_rois == 0) return 0.0;
    
    double sum_sq = 0.0;
    for (uint32_t i = 0; i < pc->n_rois; i++) {
        double sig = pc->rois[i]->sigma_line;
        if (sig > 0) {
            sum_sq += sig * sig;
        }
    }
    
    pc->significance = sqrt(sum_sq);
    
    /* Also compute combined line parameters */
    double sum_theta = 0.0, sum_rho = 0.0;
    double sum_weight = 0.0;
    
    for (uint32_t i = 0; i < pc->n_rois; i++) {
        double w = pc->rois[i]->sigma_line;
        if (w > 0) {
            sum_theta += w * pc->rois[i]->theta_line;
            sum_rho += w * pc->rois[i]->rho_line;
            sum_weight += w;
        }
    }
    
    if (sum_weight > 0) {
        pc->theta_line = sum_theta / sum_weight;
        pc->rho_line = sum_rho / sum_weight;
    }
    
    /* Compute proper motion statistics */
    if (pc->n_line_stars > 0) {
        double sum_mu_phi = 0.0, sum_mu_lambda = 0.0;
        for (uint32_t i = 0; i < pc->n_line_stars; i++) {
            sum_mu_phi += pc->line_stars[i]->mu_phi;
            sum_mu_lambda += pc->line_stars[i]->mu_lambda;
        }
        pc->mean_mu_phi = sum_mu_phi / pc->n_line_stars;
        pc->mean_mu_lambda = sum_mu_lambda / pc->n_line_stars;
        
        double var_mu_phi = 0.0, var_mu_lambda = 0.0;
        for (uint32_t i = 0; i < pc->n_line_stars; i++) {
            var_mu_phi += pow(pc->line_stars[i]->mu_phi - pc->mean_mu_phi, 2);
            var_mu_lambda += pow(pc->line_stars[i]->mu_lambda - pc->mean_mu_lambda, 2);
        }
        pc->std_mu_phi = sqrt(var_mu_phi / pc->n_line_stars);
        pc->std_mu_lambda = sqrt(var_mu_lambda / pc->n_line_stars);
    }
    
    /* Check if edge cluster */
    pc->is_edge = (fabs(pc->rho_line) > EDGE_RADIUS_DEG);
    
    return pc->significance;
}

/* ============================================================================
 * Protostream Management
 * ============================================================================ */

/**
 * @brief Create a new protostream
 * 
 * @return Newly allocated Protostream, or NULL on error
 */
Protostream *protostream_create(void) {
    Protostream *ps = (Protostream *)calloc(1, sizeof(Protostream));
    if (!ps) return NULL;
    
    ps->clusters = (Protocluster **)malloc(MAX_PROTOCLUSTERS * sizeof(Protocluster *));
    if (!ps->clusters) {
        free(ps);
        return NULL;
    }
    
    ps->n_clusters = 0;
    ps->passes_cuts = true;
    
    return ps;
}

/**
 * @brief Destroy a protostream
 * 
 * @param ps Protostream to destroy
 */
void protostream_destroy(Protostream *ps) {
    if (ps) {
        free(ps->clusters);
        free(ps);
    }
}

/**
 * @brief Add a protocluster to a protostream
 * 
 * @param ps Protostream
 * @param pc Protocluster to add
 * @return 0 on success, -1 on error
 */
int protostream_add_cluster(Protostream *ps, Protocluster *pc) {
    if (!ps || !pc) return -1;
    
    ps->clusters[ps->n_clusters++] = pc;
    
    /* Update significance (highest among clusters) */
    if (pc->significance > ps->significance) {
        ps->significance = pc->significance;
    }
    
    /* Update edge status */
    if (pc->is_edge) {
        ps->is_edge = true;
    }
    
    return 0;
}

/* Isochrone-based f_dim cut removed - algorithm uses only positions and proper motions */

/* ============================================================================
 * Stream Candidate Management
 * ============================================================================ */

/**
 * @brief Create a new stream candidate
 * 
 * @return Newly allocated StreamCandidate, or NULL on error
 */
StreamCandidate *stream_candidate_create(void) {
    StreamCandidate *sc = (StreamCandidate *)calloc(1, sizeof(StreamCandidate));
    if (!sc) return NULL;
    
    sc->protostreams = (Protostream **)malloc(MAX_PROTOSTREAMS * sizeof(Protostream *));
    sc->all_stars = (Star **)malloc(100000 * sizeof(Star *));
    
    if (!sc->protostreams || !sc->all_stars) {
        free(sc->protostreams);
        free(sc->all_stars);
        free(sc);
        return NULL;
    }
    
    sc->n_protostreams = 0;
    sc->n_stars = 0;
    memset(sc->name, 0, sizeof(sc->name));
    
    return sc;
}

/**
 * @brief Destroy a stream candidate
 * 
 * @param sc StreamCandidate to destroy
 */
void stream_candidate_destroy(StreamCandidate *sc) {
    if (sc) {
        free(sc->protostreams);
        free(sc->all_stars);
        free(sc);
    }
}

/**
 * @brief Add a protostream to a stream candidate
 * 
 * @param sc StreamCandidate
 * @param ps Protostream to add
 * @return 0 on success, -1 on error
 */
int stream_candidate_add_protostream(StreamCandidate *sc, Protostream *ps) {
    if (!sc || !ps) return -1;
    
    /* Bounds check for protostreams array */
    if (sc->n_protostreams >= MAX_PROTOSTREAMS) return -1;
    
    sc->protostreams[sc->n_protostreams++] = ps;
    
    /* Add unique stars with bounds checking */
    #define MAX_CANDIDATE_STARS 100000
    for (uint32_t c = 0; c < ps->n_clusters; c++) {
        Protocluster *pc = ps->clusters[c];
        if (!pc) continue;
        
        for (uint32_t i = 0; i < pc->n_line_stars; i++) {
            if (!pc->line_stars[i]) continue;
            if (sc->n_stars >= MAX_CANDIDATE_STARS) break;
            
            bool found = false;
            for (uint32_t j = 0; j < sc->n_stars; j++) {
                if (sc->all_stars[j] == pc->line_stars[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                sc->all_stars[sc->n_stars++] = pc->line_stars[i];
            }
        }
    }
    #undef MAX_CANDIDATE_STARS
    
    return 0;
}

/**
 * @brief Compute the final significance of a stream candidate
 * 
 * Uses quadrature sum of protostream significances.
 * 
 * @param sc StreamCandidate
 * @return Combined significance
 */
double stream_candidate_compute_significance(StreamCandidate *sc) {
    if (!sc || sc->n_protostreams == 0) return 0.0;
    
    double sum_sq = 0.0;
    double sum_ra = 0.0, sum_dec = 0.0;
    double sum_pmra = 0.0, sum_pmdec = 0.0;
    
    for (uint32_t i = 0; i < sc->n_protostreams; i++) {
        if (!sc->protostreams[i]) continue;
        double sig = sc->protostreams[i]->significance;
        sum_sq += sig * sig;
    }
    
    sc->significance = sqrt(sum_sq);
    
    /* Compute position summary */
    if (sc->n_stars > 0 && sc->all_stars) {
        uint32_t valid_stars = 0;
        for (uint32_t i = 0; i < sc->n_stars; i++) {
            if (!sc->all_stars[i]) continue;
            sum_ra += sc->all_stars[i]->ra;
            sum_dec += sc->all_stars[i]->dec;
            sum_pmra += sc->all_stars[i]->pmra;
            sum_pmdec += sc->all_stars[i]->pmdec;
            valid_stars++;
        }
        
        if (valid_stars > 0) {
            sc->mean_ra = sum_ra / valid_stars;
            sc->mean_dec = sum_dec / valid_stars;
            sc->mean_pmra = sum_pmra / valid_stars;
            sc->mean_pmdec = sum_pmdec / valid_stars;
            
            /* Convert to galactic */
            double l, b;
            extern void equatorial_to_galactic(double ra, double dec, double *l, double *b);
            equatorial_to_galactic(sc->mean_ra, sc->mean_dec, &l, &b);
            sc->mean_l = l;
            sc->mean_b = b;
        }
    }
    
    return sc->significance;
}

/* ============================================================================
 * Clustering Algorithms
 * ============================================================================ */

/**
 * @brief Check if two ROIs can be clustered together
 * 
 * ROIs must be from different (independent) signal regions and have
 * compatible proper motions.
 * 
 * @param roi1 First ROI
 * @param roi2 Second ROI
 * @param pm_threshold Maximum proper motion distance (mas/yr)
 * @return true if ROIs can be clustered
 */
static bool rois_can_cluster(ROI *roi1, ROI *roi2, double pm_threshold) {
    /* ROIs from same SR cannot be clustered (not independent) */
    if (roi1->sr_id == roi2->sr_id) {
        return false;
    }
    
    /* Check proper motion compatibility */
    /* Using line stars' proper motions */
    double mu_phi_1 = 0.0, mu_lambda_1 = 0.0;
    double mu_phi_2 = 0.0, mu_lambda_2 = 0.0;
    
    for (uint32_t i = 0; i < roi1->n_selected && i < 10; i++) {
        mu_phi_1 += roi1->selected_stars[i]->mu_phi;
        mu_lambda_1 += roi1->selected_stars[i]->mu_lambda;
    }
    int n1 = (roi1->n_selected < 10) ? roi1->n_selected : 10;
    mu_phi_1 /= n1;
    mu_lambda_1 /= n1;
    
    for (uint32_t i = 0; i < roi2->n_selected && i < 10; i++) {
        mu_phi_2 += roi2->selected_stars[i]->mu_phi;
        mu_lambda_2 += roi2->selected_stars[i]->mu_lambda;
    }
    int n2 = (roi2->n_selected < 10) ? roi2->n_selected : 10;
    mu_phi_2 /= n2;
    mu_lambda_2 /= n2;
    
    double pm_dist = sqrt(pow(mu_phi_1 - mu_phi_2, 2) + 
                          pow(mu_lambda_1 - mu_lambda_2, 2));
    
    return (pm_dist < pm_threshold);
}

/**
 * @brief Cluster ROIs into protoclusters
 * 
 * Uses a greedy clustering algorithm that combines ROIs with compatible
 * lines and proper motions.
 * 
 * @param rois Array of ROIs
 * @param n_rois Number of ROIs
 * @param clusters Output array of protoclusters
 * @param n_clusters Output number of protoclusters
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int cluster_rois_to_protoclusters(ROI **rois, int n_rois,
                                  Protocluster **clusters, int *n_clusters,
                                  const Config *cfg) {
    if (!rois || !clusters || !n_clusters || !cfg) {
        return -1;
    }
    
    /* Track which ROIs have been assigned */
    bool *assigned = (bool *)calloc(n_rois, sizeof(bool));
    if (!assigned) return -1;
    
    *n_clusters = 0;
    double pm_threshold = 5.0;  /* mas/yr - increased for robustness */
    
    for (int i = 0; i < n_rois; i++) {
        if (assigned[i]) continue;
        if (rois[i]->sigma_line < 2.0) continue;  /* Lower threshold for testing */
        
        /* Start new cluster with this ROI */
        Protocluster *pc = protocluster_create();
        if (!pc) {
            free(assigned);
            return -1;
        }
        
        protocluster_add_roi(pc, rois[i]);
        pc->patch_id = rois[i]->patch_id;
        assigned[i] = true;
        
        /* Try to add more ROIs */
        bool added = true;
        while (added) {
            added = false;
            for (int j = 0; j < n_rois; j++) {
                if (assigned[j]) continue;
                
                /* Check if this ROI can be added */
                bool can_add = false;
                for (uint32_t k = 0; k < pc->n_rois; k++) {
                    if (rois_can_cluster(pc->rois[k], rois[j], pm_threshold)) {
                        /* Check line compatibility */
                        if (lines_compatible(pc->rois[k]->theta_line, pc->rois[k]->rho_line,
                                           rois[j]->theta_line, rois[j]->rho_line,
                                           cfg->theta_merge_threshold, cfg->rho_merge_threshold)) {
                            can_add = true;
                            break;
                        }
                    }
                }
                
                if (can_add) {
                    protocluster_add_roi(pc, rois[j]);
                    assigned[j] = true;
                    added = true;
                }
            }
        }
        
        /* Compute significance */
        protocluster_compute_significance(pc);
        
        /* Keep if above threshold */
        if (pc->significance >= cfg->protocluster_sig_cut) {
            pc->cluster_id = *n_clusters;
            clusters[(*n_clusters)++] = pc;
        } else {
            protocluster_destroy(pc);
        }
    }
    
    free(assigned);
    return 0;
}

/**
 * @brief Deduplicate protoclusters into protostreams
 * 
 * Protoclusters with >40% star overlap are considered duplicates
 * and grouped into a protostream.
 * 
 * @param clusters Array of protoclusters
 * @param n_clusters Number of protoclusters
 * @param protostreams Output array of protostreams
 * @param n_protostreams Output number of protostreams
 * @return 0 on success, error code otherwise
 */
int deduplicate_protoclusters(Protocluster **clusters, int n_clusters,
                             Protostream **protostreams, int *n_protostreams) {
    if (!clusters || !protostreams || !n_protostreams) {
        return -1;
    }
    
    /* Build overlap matrix */
    double *overlap = (double *)calloc(n_clusters * n_clusters, sizeof(double));
    if (!overlap) return -1;
    
    for (int i = 0; i < n_clusters; i++) {
        if (!clusters[i]) continue;
        for (int j = i + 1; j < n_clusters; j++) {
            if (!clusters[j]) continue;
            if (clusters[i]->n_line_stars == 0 || clusters[j]->n_line_stars == 0) continue;
            
            double ov = star_set_overlap(clusters[i]->line_stars, clusters[i]->n_line_stars,
                                        clusters[j]->line_stars, clusters[j]->n_line_stars);
            overlap[i * n_clusters + j] = ov;
            overlap[j * n_clusters + i] = ov;
        }
    }
    
    /* Union-Find for grouping */
    int *parent = (int *)malloc(n_clusters * sizeof(int));
    if (!parent) {
        free(overlap);
        return -1;
    }
    for (int i = 0; i < n_clusters; i++) {
        parent[i] = i;
    }
    
    /* Group overlapping clusters */
    for (int i = 0; i < n_clusters; i++) {
        for (int j = i + 1; j < n_clusters; j++) {
            if (overlap[i * n_clusters + j] > DUPLICATE_OVERLAP) {
                uf_unite(parent, i, j);
            }
        }
    }
    
    /* Create protostreams from groups */
    *n_protostreams = 0;
    bool *processed = (bool *)calloc(n_clusters, sizeof(bool));
    if (!processed) {
        free(overlap);
        free(parent);
        return -1;
    }
    
    for (int i = 0; i < n_clusters; i++) {
        int root = uf_find(parent, i);
        if (processed[root]) continue;
        
        Protostream *ps = protostream_create();
        if (!ps) {
            free(overlap);
            free(parent);
            free(processed);
            return -1;
        }
        
        /* Add all clusters in this group */
        for (int j = 0; j < n_clusters; j++) {
            if (uf_find(parent, j) == root) {
                protostream_add_cluster(ps, clusters[j]);
            }
        }
        
        ps->protostream_id = *n_protostreams;
        ps->patch_id = clusters[i]->patch_id;
        processed[root] = true;
        
        protostreams[(*n_protostreams)++] = ps;
    }
    
    free(overlap);
    free(parent);
    free(processed);
    
    return 0;
}

/**
 * @brief Merge protostreams across patches into stream candidates
 * 
 * Protostreams from adjacent/overlapping patches are merged if they
 * share stars and have compatible line parameters.
 * 
 * @param protostreams Array of all protostreams from all patches
 * @param n_protostreams Number of protostreams
 * @param candidates Output array of stream candidates
 * @param n_candidates Output number of candidates
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int merge_protostreams_across_patches(Protostream **protostreams, int n_protostreams,
                                      StreamCandidate **candidates, int *n_candidates,
                                      const Config *cfg) {
    if (!protostreams || !candidates || !n_candidates || !cfg) {
        return -1;
    }
    
    /* Union-Find for grouping */
    int *parent = (int *)malloc(n_protostreams * sizeof(int));
    if (!parent) {
        return -1;
    }
    for (int i = 0; i < n_protostreams; i++) {
        parent[i] = i;
    }
    
    printf("  Checking %d protostreams for merging...\n", n_protostreams);
    
    /* Check each pair of protostreams */
    for (int i = 0; i < n_protostreams; i++) {
        /* Null check */
        if (!protostreams[i]) continue;
        
        for (int j = i + 1; j < n_protostreams; j++) {
            /* Null check */
            if (!protostreams[j]) continue;
            
            /* Must be from different patches */
            if (protostreams[i]->patch_id == protostreams[j]->patch_id) {
                continue;
            }
            
            /* Check if clusters overlap and have compatible lines */
            bool should_merge = false;
            
            for (uint32_t ci = 0; ci < protostreams[i]->n_clusters && !should_merge; ci++) {
                Protocluster *pci = protostreams[i]->clusters[ci];
                if (!pci || pci->n_line_stars == 0) continue;
                
                for (uint32_t cj = 0; cj < protostreams[j]->n_clusters && !should_merge; cj++) {
                    Protocluster *pcj = protostreams[j]->clusters[cj];
                    if (!pcj || pcj->n_line_stars == 0) continue;
                    
                    /* Check star overlap */
                    double ov = star_set_overlap(pci->line_stars, pci->n_line_stars,
                                                pcj->line_stars, pcj->n_line_stars);
                    
                    if (ov > cfg->line_overlap_threshold) {
                        /* Check line compatibility */
                        if (lines_compatible(pci->theta_line, pci->rho_line,
                                           pcj->theta_line, pcj->rho_line,
                                           cfg->theta_merge_threshold, cfg->rho_merge_threshold)) {
                            should_merge = true;
                        }
                    }
                }
            }
            
            if (should_merge) {
                uf_unite(parent, i, j);
            }
        }
        
        /* Progress indicator for large datasets */
        if (i > 0 && i % 100 == 0) {
            printf("    Processed %d/%d protostream pairs...\n", i, n_protostreams);
        }
    }
    
    printf("  Creating stream candidates from merged groups...\n");
    
    /* Create stream candidates from groups */
    *n_candidates = 0;
    bool *processed = (bool *)calloc(n_protostreams, sizeof(bool));
    if (!processed) {
        free(parent);
        return -1;
    }
    
    for (int i = 0; i < n_protostreams; i++) {
        int root = uf_find(parent, i);
        if (processed[root]) continue;
        
        StreamCandidate *sc = stream_candidate_create();
        if (!sc) {
            free(parent);
            free(processed);
            return -1;
        }
        
        /* Add all protostreams in this group */
        bool has_non_edge = false;
        for (int j = 0; j < n_protostreams; j++) {
            if (uf_find(parent, j) == root) {
                stream_candidate_add_protostream(sc, protostreams[j]);
                if (!protostreams[j]->is_edge) {
                    has_non_edge = true;
                }
            }
        }
        
        /* Remove singleton edge protostreams */
        if (sc->n_protostreams == 1 && protostreams[i]->is_edge) {
            stream_candidate_destroy(sc);
            processed[root] = true;
            continue;
        }
        
        /* Mark as edge if all protostreams are edge (and no non-edge) */
        if (!has_non_edge && sc->n_protostreams > 1) {
            /* Could be suspicious, but keep for now */
        }
        
        sc->candidate_id = *n_candidates;
        stream_candidate_compute_significance(sc);
        processed[root] = true;
        
        candidates[(*n_candidates)++] = sc;
    }
    
    free(parent);
    free(processed);
    
    return 0;
}

/* Isochrone-based apply_f_dim_cut removed - algorithm uses only positions and proper motions */

