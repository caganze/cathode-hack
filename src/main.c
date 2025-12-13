/**
 * @file main.c
 * @brief Main entry point for the stream detection program
 * 
 * Via Machinae-inspired stellar stream detection for DESI DR1 MWS data
 * 
 * Based on arXiv:2509.08064v1
 * "Via Machinae 3.0: A search for stellar streams in Gaia with the CATHODE algorithm"
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>  /* For strcasecmp */
#include <getopt.h>
#include <time.h>
#include <math.h>     /* For cos, fabs */
#include <sys/stat.h>
#include <errno.h>
#include "stream_detect.h"

/* Version string */
#define VERSION "1.0.0"

/* ============================================================================
 * Configuration Management
 * ============================================================================ */

/**
 * @brief Create default configuration
 */
Config *config_create_default(void) {
    Config *cfg = (Config *)calloc(1, sizeof(Config));
    if (!cfg) return NULL;
    
    /* Patch parameters */
    cfg->patch_radius = PATCH_RADIUS_DEG;
    cfg->min_galactic_b = 20.0;  /* Exclude |b| < 20 degrees */
    
    /* Signal region parameters */
    cfg->sr_width = SR_WIDTH_MAS_YR;
    cfg->sr_step = SR_STEP_MAS_YR;
    cfg->min_stars_sr = MIN_STARS_IN_SR;
    cfg->max_stars_sr = MAX_STARS_IN_SR;
    
    /* ROI parameters */
    cfg->roi_width = ROI_WIDTH_MAS_YR;
    cfg->top_n_anomalous = TOP_N_ANOMALOUS;
    
    /* Fiducial cuts - position-based only (no isochrone) */
    cfg->dist_limit = FIDUCIAL_DIST_DEG;
    cfg->z_min_kpc = 2.0;
    
    /* Significance cuts - aggressive for dense DESI backgrounds */
    cfg->protocluster_sig_cut = 8.0;  /* Higher threshold for real data with dense background */
    cfg->edge_radius = EDGE_RADIUS_DEG;
    
    /* Clustering parameters */
    cfg->duplicate_overlap = DUPLICATE_OVERLAP;
    cfg->line_overlap_threshold = 0.1;
    cfg->theta_merge_threshold = 14.0 * (3.14159265 / 180.0);  /* 14 degrees in radians */
    cfg->rho_merge_threshold = 1.0;  /* degrees */
    
    /* Sky region selection (default: full sky) */
    cfg->ra_min = 0.0;
    cfg->ra_max = 360.0;
    cfg->dec_min = -90.0;
    cfg->dec_max = 90.0;
    
    /* Distance cuts (default: no cut) */
    cfg->dist_min_kpc = 0.0;
    cfg->dist_max_kpc = 1000.0;
    
    /* Output */
    strcpy(cfg->output_dir, "./output");
    cfg->verbose = true;
    
    printf("Configuration:\n");
    printf("  patch_radius: %.1f deg\n", cfg->patch_radius);
    printf("  min_galactic_b: %.1f deg\n", cfg->min_galactic_b);
    printf("  sr_width: %.1f mas/yr\n", cfg->sr_width);
    printf("  protocluster_sig_cut: %.1f\n", cfg->protocluster_sig_cut);
    printf("  top_n_anomalous: %d\n", cfg->top_n_anomalous);
    
    return cfg;
}

/**
 * @brief Destroy configuration
 */
void config_destroy(Config *cfg) {
    free(cfg);
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Create directory if it doesn't exist
 */
static int ensure_directory(const char *path) {
    struct stat st;
    if (stat(path, &st) == 0) {
        return 0;  /* Already exists */
    }
    
    if (mkdir(path, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "Error: Cannot create directory %s: %s\n", 
                path, strerror(errno));
        return -1;
    }
    
    return 0;
}

/**
 * @brief Print usage information
 */
static void print_usage(const char *progname) {
    printf("Usage: %s [options] <input_file>\n\n", progname);
    printf("Via Machinae-inspired stellar stream detection for DESI DR1 MWS data\n\n");
    printf("Options:\n");
    printf("  -h, --help              Show this help message\n");
    printf("  -v, --version           Show version\n");
    printf("  -o, --output DIR        Output directory (default: ./output)\n");
    printf("  -p, --patch-radius R    Patch radius in degrees (default: 1.78)\n");
    printf("  -b, --min-b B           Minimum galactic |b| (default: 20.0)\n");
    printf("  -s, --sig-cut S         Protocluster significance cut (default: 8.0)\n");
    printf("  -k, --knn K             Use KNN density with K neighbors (default: KDE)\n");
    printf("  --histogram             Use histogram density (faster, less accurate)\n");
    printf("  -q, --quiet             Quiet mode (less output)\n\n");
    printf("Sky region selection:\n");
    printf("  --ra-min RA             Minimum RA in degrees (default: 0)\n");
    printf("  --ra-max RA             Maximum RA in degrees (default: 360)\n");
    printf("  --dec-min DEC           Minimum Dec in degrees (default: -90)\n");
    printf("  --dec-max DEC           Maximum Dec in degrees (default: 90)\n\n");
    printf("Distance cuts:\n");
    printf("  --dist-min D            Minimum distance in kpc (default: 0)\n");
    printf("  --dist-max D            Maximum distance in kpc (default: 1000)\n\n");
    printf("Examples:\n");
    printf("  %s data.csv --ra-min 100 --ra-max 200 --dec-min -30 --dec-max 30\n", progname);
    printf("  %s data.csv --dist-min 5 --dist-max 50 --sig-cut 10\n\n", progname);
    printf("Input file can be CSV or FITS format.\n");
    printf("For FITS support, compile with -DUSE_CFITSIO and link with -lcfitsio\n");
}

/* ============================================================================
 * Main Pipeline
 * ============================================================================ */

/**
 * @brief Run the full stream detection pipeline
 */
int run_stream_detection(Star *stars, uint32_t n_stars, const Config *cfg,
                         StreamCandidate ***candidates, int *n_candidates) {
    if (!stars || n_stars == 0 || !cfg || !candidates || !n_candidates) {
        return -1;
    }
    
    time_t start_time = time(NULL);
    
    printf("\n========================================\n");
    printf("Stream Detection Pipeline\n");
    printf("========================================\n");
    printf("Input: %u stars\n", n_stars);
    printf("Patch radius: %.1f degrees\n", cfg->patch_radius);
    printf("Min galactic |b|: %.1f degrees\n", cfg->min_galactic_b);
    printf("Protocluster significance cut: %.2f\n", cfg->protocluster_sig_cut);
    printf("\n");
    
    /* ========================================================================
     * Step 1: Generate patches and assign stars
     * ======================================================================== */
    printf("Step 1: Generating sky patches...\n");
    
    Patch **patches = NULL;
    int n_patches = 0;
    
    /* Allocate patch pointer array */
    patches = (Patch **)malloc(MAX_PATCHES * sizeof(Patch *));
    if (!patches) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return -1;
    }
    
    /* Generate patches manually */
    /* 10 sq.deg patches with some overlap */
    n_patches = 0;
    double spacing = cfg->patch_radius * 1.5;  /* Some overlap between patches */
    
    for (double dec = -80; dec <= 80; dec += spacing) {
        double cos_dec = cos(dec * 3.14159265 / 180.0);
        if (cos_dec < 0.1) cos_dec = 0.1;
        double ra_spacing = spacing / cos_dec;
        
        for (double ra = 0; ra < 360; ra += ra_spacing) {
            double l, b;
            extern void equatorial_to_galactic(double ra, double dec, double *l, double *b);
            equatorial_to_galactic(ra, dec, &l, &b);
            
            if (fabs(b) < cfg->min_galactic_b) continue;
            
            /* Skip LMC/SMC */
            extern double angular_separation(double ra1, double dec1, double ra2, double dec2);
            if (angular_separation(ra, dec, 80.0, -69.0) < 10.0) continue;
            if (angular_separation(ra, dec, 13.0, -73.0) < 5.0) continue;
            
            if (n_patches >= MAX_PATCHES) break;
            
            patches[n_patches] = patch_create(ra, dec, cfg->patch_radius);
            if (patches[n_patches]) {
                patches[n_patches]->patch_id = n_patches;
                n_patches++;
            }
        }
        if (n_patches >= MAX_PATCHES) break;
    }
    
    printf("  Created %d patches\n", n_patches);
    
    /* Assign stars to patches */
    printf("  Assigning stars to patches...\n");
    
    int total_assigned = 0;
    for (uint32_t i = 0; i < n_stars; i++) {
        Star *star = &stars[i];
        
        for (int p = 0; p < n_patches; p++) {
            double dist = angular_separation(star->ra, star->dec,
                                           patches[p]->center_ra, 
                                           patches[p]->center_dec);
            
            if (dist <= patches[p]->radius) {
                extern int patch_add_star(Patch *patch, const Star *star);
                patch_add_star(patches[p], star);
                total_assigned++;
            }
        }
        
        if (i % 500000 == 0 && i > 0) {
            printf("    Processed %u stars...\n", i);
        }
    }
    
    printf("  Total star-patch assignments: %d\n\n", total_assigned);
    
    /* ========================================================================
     * Step 2: Process each patch
     * ======================================================================== */
    printf("Step 2: Processing patches...\n");
    
    /* Collect all protostreams */
    Protostream **all_protostreams = (Protostream **)malloc(
        MAX_PROTOSTREAMS * n_patches * sizeof(Protostream *));
    int total_protostreams = 0;
    
    for (int p = 0; p < n_patches; p++) {
        Patch *patch = patches[p];
        
        if (patch->n_stars < 50) {
            continue;  /* Skip patches with too few stars */
        }
        
        if (cfg->verbose && p % 10 == 0) {
            printf("  Processing patch %d/%d (stars: %u)\n", 
                   p + 1, n_patches, patch->n_stars);
        }
        
        /* Apply fiducial cuts - but be lenient for testing */
        extern int apply_fiducial_cuts(Patch *patch, const Config *cfg);
        int n_passing = apply_fiducial_cuts(patch, cfg);
        
        if (cfg->verbose && p % 50 == 0) {
            printf("    Patch %d: %d/%u stars pass fiducial cuts\n", 
                   p, n_passing, patch->n_stars);
        }
        
        /* Define signal regions and compute anomaly scores */
        ROI **patch_rois = (ROI **)malloc(MAX_ROIS_PER_PATCH * sizeof(ROI *));
        int n_patch_rois = 0;
        int sr_id = 0;
        int n_srs_created = 0;
        int n_rois_created = 0;
        
        /* Iterate over proper motion windows */
        /* PM range covers typical halo/stream motions */
        double mu_min = -30.0;
        double mu_max = 30.0;
        
        int n_high_sig_rois = 0;
        double max_roi_sig = 0.0;
        
        /* Scan PM space with 0.5 mas/yr windows */
        for (double mu = mu_min; mu < mu_max - cfg->sr_width; mu += cfg->sr_step) {
            /* Try both pmdec and pmra as the "resonant" feature */
            for (int use_lambda = 0; use_lambda <= 1; use_lambda++) {
                SignalRegion *sr = sr_create(patch, mu, mu + cfg->sr_width, use_lambda);
                if (!sr) continue;
                
                /* Check star counts - adjusted for smaller patches */
                if (sr->n_sr_stars < (uint32_t)cfg->min_stars_sr || 
                    sr->n_sr_stars > (uint32_t)cfg->max_stars_sr) {
                    sr_destroy(sr);
                    continue;
                }
                
                n_srs_created++;
                
                /* Compute anomaly scores using KNN */
                int k = (sr->n_sb_stars < 15) ? sr->n_sb_stars : 15;
                if (k < 3) k = 3;
                compute_anomaly_scores_knn(patch, sr, k, cfg);
                
                /* Create ROIs within this SR */
                for (double mu_p = mu_min; mu_p < mu_max - cfg->roi_width; mu_p += cfg->sr_step) {
                    ROI *roi = roi_create(sr, mu_p, mu_p + cfg->roi_width);
                    if (!roi) continue;
                    
                    n_rois_created++;
                    
                    /* Find line */
                    find_line_in_roi(roi, cfg);
                    
                    /* Track max significance */
                    if (roi->sigma_line > max_roi_sig) {
                        max_roi_sig = roi->sigma_line;
                    }
                    
                    /* Keep if significant - aggressive cutoff for dense backgrounds */
                    if (roi->sigma_line >= 4.0) {
                        roi->sr_id = sr_id;
                        roi->patch_id = p;
                        n_high_sig_rois++;
                        
                        if (n_patch_rois < MAX_ROIS_PER_PATCH) {
                            patch_rois[n_patch_rois++] = roi;
                        } else {
                            roi_destroy(roi);
                        }
                    } else {
                        roi_destroy(roi);
                    }
                }
                
                sr_id++;
                sr_destroy(sr);
            }
        }
        
        /* Debug output for patches with potential detections */
        if (cfg->verbose && (n_high_sig_rois > 0 || p % 50 == 0)) {
            printf("    Patch %d: %d SRs, %d ROIs, %d sig>=3.0, max_sig=%.2f\n",
                   p, n_srs_created, n_rois_created, n_high_sig_rois, max_roi_sig);
        }
        
        if (cfg->verbose && p % 50 == 0) {
            printf("    Patch %d: %d SRs, %d ROIs created, %d significant ROIs\n",
                   p, n_srs_created, n_rois_created, n_patch_rois);
        }
        
        /* Cluster ROIs into protoclusters */
        if (n_patch_rois > 0) {
            if (cfg->verbose && p % 50 == 0) {
                printf("    Patch %d: Clustering %d significant ROIs\n", p, n_patch_rois);
            }
            
            Protocluster **protoclusters = (Protocluster **)malloc(
                MAX_PROTOCLUSTERS * sizeof(Protocluster *));
            int n_protoclusters = 0;
            
            cluster_rois_to_protoclusters(patch_rois, n_patch_rois,
                                         protoclusters, &n_protoclusters, cfg);
            
            if (cfg->verbose && n_protoclusters > 0) {
                printf("    Patch %d: Found %d protoclusters\n", p, n_protoclusters);
            }
            
            /* Deduplicate into protostreams */
            if (n_protoclusters > 0) {
                Protostream **patch_protostreams = (Protostream **)malloc(
                    MAX_PROTOSTREAMS * sizeof(Protostream *));
                int n_patch_protostreams = 0;
                
                deduplicate_protoclusters(protoclusters, n_protoclusters,
                                         patch_protostreams, &n_patch_protostreams);
                
                /* Collect valid protostreams (no isochrone cut) */
                for (int ps = 0; ps < n_patch_protostreams; ps++) {
                    /* Use aggressive significance cut for dense backgrounds */
                    if (patch_protostreams[ps]->significance >= cfg->protocluster_sig_cut) {
                        patch_protostreams[ps]->passes_cuts = true;
                        all_protostreams[total_protostreams++] = patch_protostreams[ps];
                        
                        if (cfg->verbose) {
                            printf("    Patch %d: Protostream with sig=%.2f\n",
                                   p, patch_protostreams[ps]->significance);
                        }
                    } else {
                        patch_protostreams[ps]->passes_cuts = false;
                        protostream_destroy(patch_protostreams[ps]);
                    }
                }
                
                free(patch_protostreams);
            }
            
            free(protoclusters);
        }
        
        /* Cleanup ROIs */
        for (int r = 0; r < n_patch_rois; r++) {
            /* ROIs are referenced by protoclusters, handle carefully */
        }
        free(patch_rois);
    }
    
    printf("  Found %d protostreams passing cuts\n\n", total_protostreams);
    
    /* ========================================================================
     * Step 3: Merge protostreams across patches
     * ======================================================================== */
    printf("Step 3: Merging protostreams across patches...\n");
    
    *candidates = (StreamCandidate **)malloc(MAX_STREAM_CANDIDATES * sizeof(StreamCandidate *));
    *n_candidates = 0;
    
    if (total_protostreams > 0) {
        merge_protostreams_across_patches(all_protostreams, total_protostreams,
                                         *candidates, n_candidates, cfg);
    }
    
    /* Compute final significances and sort */
    for (int i = 0; i < *n_candidates; i++) {
        if ((*candidates)[i]) {
            stream_candidate_compute_significance((*candidates)[i]);
        }
    }
    
    /* Sort by significance (with null checks) */
    for (int i = 0; i < *n_candidates - 1; i++) {
        for (int j = i + 1; j < *n_candidates; j++) {
            if (!(*candidates)[i] || !(*candidates)[j]) continue;
            if ((*candidates)[j]->significance > (*candidates)[i]->significance) {
                StreamCandidate *tmp = (*candidates)[i];
                (*candidates)[i] = (*candidates)[j];
                (*candidates)[j] = tmp;
            }
        }
    }
    
    printf("  Found %d stream candidates\n\n", *n_candidates);
    
    time_t end_time = time(NULL);
    printf("Pipeline completed in %ld seconds\n", end_time - start_time);
    printf("========================================\n\n");
    
    /* Note: We don't destroy patches here because stream candidates
     * still reference stars within patches. Memory will be freed on exit.
     * For long-running applications, proper reference counting would be needed.
     */
    free(patches);
    free(all_protostreams);
    
    return 0;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

int main(int argc, char *argv[]) {
    /* Default options */
    Config *cfg = config_create_default();
    char input_file[512] = "";
    int use_knn = 0;
    int knn_k = 20;
    int use_histogram = 0;
    
    /* Parse command line options */
    static struct option long_options[] = {
        {"help",         no_argument,       0, 'h'},
        {"version",      no_argument,       0, 'v'},
        {"output",       required_argument, 0, 'o'},
        {"patch-radius", required_argument, 0, 'p'},
        {"min-b",        required_argument, 0, 'b'},
        {"sig-cut",      required_argument, 0, 's'},
        {"knn",          required_argument, 0, 'k'},
        {"histogram",    no_argument,       0, 'H'},
        {"quiet",        no_argument,       0, 'q'},
        {"ra-min",       required_argument, 0, 1001},
        {"ra-max",       required_argument, 0, 1002},
        {"dec-min",      required_argument, 0, 1003},
        {"dec-max",      required_argument, 0, 1004},
        {"dist-min",     required_argument, 0, 1005},
        {"dist-max",     required_argument, 0, 1006},
        {0, 0, 0, 0}
    };
    
    int opt;
    while ((opt = getopt_long(argc, argv, "hvo:p:b:s:k:Hq", long_options, NULL)) != -1) {
        switch (opt) {
            case 'h':
                print_usage(argv[0]);
                config_destroy(cfg);
                return 0;
            case 'v':
                printf("stream_detect version %s\n", VERSION);
                config_destroy(cfg);
                return 0;
            case 'o':
                strncpy(cfg->output_dir, optarg, sizeof(cfg->output_dir) - 1);
                break;
            case 'p':
                cfg->patch_radius = atof(optarg);
                break;
            case 'b':
                cfg->min_galactic_b = atof(optarg);
                break;
            case 's':
                cfg->protocluster_sig_cut = atof(optarg);
                break;
            case 'k':
                use_knn = 1;
                knn_k = atoi(optarg);
                break;
            case 'H':
                use_histogram = 1;
                break;
            case 'q':
                cfg->verbose = false;
                break;
            case 1001:  /* --ra-min */
                cfg->ra_min = atof(optarg);
                break;
            case 1002:  /* --ra-max */
                cfg->ra_max = atof(optarg);
                break;
            case 1003:  /* --dec-min */
                cfg->dec_min = atof(optarg);
                break;
            case 1004:  /* --dec-max */
                cfg->dec_max = atof(optarg);
                break;
            case 1005:  /* --dist-min */
                cfg->dist_min_kpc = atof(optarg);
                break;
            case 1006:  /* --dist-max */
                cfg->dist_max_kpc = atof(optarg);
                break;
            default:
                print_usage(argv[0]);
                config_destroy(cfg);
                return 1;
        }
    }
    
    /* Get input file */
    if (optind >= argc) {
        fprintf(stderr, "Error: No input file specified\n\n");
        print_usage(argv[0]);
        config_destroy(cfg);
        return 1;
    }
    strncpy(input_file, argv[optind], sizeof(input_file) - 1);
    
    printf("\n");
    printf("========================================\n");
    printf("Via Machinae Stream Detector v%s\n", VERSION);
    printf("========================================\n");
    printf("Based on arXiv:2509.08064v1\n");
    printf("Input: %s\n", input_file);
    printf("Output: %s\n", cfg->output_dir);
    printf("\n");
    
    /* Create output directory */
    if (ensure_directory(cfg->output_dir) != 0) {
        config_destroy(cfg);
        return 1;
    }
    
    /* Read input data */
    Star *stars = NULL;
    uint32_t n_stars = 0;
    int read_result;
    
    /* Check file extension */
    char *ext = strrchr(input_file, '.');
    if (ext && (strcasecmp(ext, ".fits") == 0 || strcasecmp(ext, ".fit") == 0)) {
        read_result = read_desi_mws_fits(input_file, &stars, &n_stars);
    } else {
        read_result = read_stars_csv(input_file, &stars, &n_stars);
    }
    
    if (read_result != 0 || n_stars == 0) {
        fprintf(stderr, "Error: Failed to read input file\n");
        config_destroy(cfg);
        return 1;
    }
    
    /* Apply RA/Dec/distance filters */
    bool has_sky_cut = (cfg->ra_min > 0.0 || cfg->ra_max < 360.0 || 
                        cfg->dec_min > -90.0 || cfg->dec_max < 90.0);
    bool has_dist_cut = (cfg->dist_min_kpc > 0.0 || cfg->dist_max_kpc < 999.0);
    
    if (has_sky_cut || has_dist_cut) {
        printf("\nApplying selection cuts:\n");
        if (has_sky_cut) {
            printf("  RA range:  [%.2f, %.2f] deg\n", cfg->ra_min, cfg->ra_max);
            printf("  Dec range: [%.2f, %.2f] deg\n", cfg->dec_min, cfg->dec_max);
        }
        if (has_dist_cut) {
            printf("  Distance:  [%.2f, %.2f] kpc\n", cfg->dist_min_kpc, cfg->dist_max_kpc);
        }
        
        uint32_t n_selected = 0;
        for (uint32_t i = 0; i < n_stars; i++) {
            Star *s = &stars[i];
            bool keep = true;
            
            /* RA cut (handle wrap-around at 360) */
            if (cfg->ra_min <= cfg->ra_max) {
                /* Normal range */
                if (s->ra < cfg->ra_min || s->ra > cfg->ra_max) keep = false;
            } else {
                /* Wrap-around (e.g., ra_min=350, ra_max=10 means 350-360 and 0-10) */
                if (s->ra < cfg->ra_min && s->ra > cfg->ra_max) keep = false;
            }
            
            /* Dec cut */
            if (s->dec < cfg->dec_min || s->dec > cfg->dec_max) keep = false;
            
            /* Distance cut */
            if (s->distance > 0) {
                if (s->distance < cfg->dist_min_kpc || s->distance > cfg->dist_max_kpc) {
                    keep = false;
                }
            }
            
            if (keep) {
                if (n_selected != i) {
                    stars[n_selected] = stars[i];
                }
                stars[n_selected].star_idx = n_selected;
                n_selected++;
            }
        }
        
        printf("  Stars before cuts: %u\n", n_stars);
        printf("  Stars after cuts:  %u\n", n_selected);
        n_stars = n_selected;
        
        if (n_stars == 0) {
            fprintf(stderr, "Error: No stars remaining after cuts\n");
            free(stars);
            config_destroy(cfg);
            return 1;
        }
    }
    
    /* Run pipeline */
    StreamCandidate **candidates = NULL;
    int n_candidates = 0;
    
    if (run_stream_detection(stars, n_stars, cfg, &candidates, &n_candidates) != 0) {
        fprintf(stderr, "Error: Pipeline failed\n");
        free(stars);
        config_destroy(cfg);
        return 1;
    }
    
    /* Write results */
    char output_path[1024];
    
    snprintf(output_path, sizeof(output_path), "%s/stream_candidates.csv", cfg->output_dir);
    write_stream_candidates_csv(output_path, candidates, n_candidates);
    
    /* Write individual stream star files for top candidates */
    int n_write = (n_candidates < 20) ? n_candidates : 20;
    for (int i = 0; i < n_write; i++) {
        if (!candidates || !candidates[i]) continue;
        if (candidates[i]->significance <= 0) continue;
        if (!candidates[i]->all_stars || candidates[i]->n_stars == 0) continue;
        
        snprintf(output_path, sizeof(output_path), "%s/stream_%03d_stars.csv", 
                 cfg->output_dir, i + 1);
        write_stream_stars_csv(output_path, candidates[i]);
    }
    
    /* Summary */
    printf("\nResults Summary:\n");
    printf("================\n");
    printf("Total stream candidates: %d\n\n", n_candidates);
    
    if (candidates && n_candidates > 0) {
        printf("Top 10 candidates by significance:\n");
        printf("%-4s %-10s %-8s %-8s %-8s %-8s\n", 
               "Rank", "Sig", "N_stars", "RA", "Dec", "l");
        printf("---- ---------- -------- -------- -------- --------\n");
        
        int n_print = (n_candidates < 10) ? n_candidates : 10;
        int rank = 0;
        for (int i = 0; i < n_print && rank < 10; i++) {
            StreamCandidate *sc = candidates[i];
            if (!sc) continue;
            if (sc->significance <= 0) continue;
            
            rank++;
            printf("%-4d %-10.2f %-8u %-8.2f %-8.2f %-8.2f\n",
                   rank, sc->significance, sc->n_stars,
                   sc->mean_ra, sc->mean_dec, sc->mean_l);
        }
    }
    
    printf("\nOutput written to: %s/\n", cfg->output_dir);
    
    /* Cleanup - note: stream candidates may reference stars in the original
     * array, so we don't deeply free them here. Memory is freed on exit.
     */
    free(candidates);
    free(stars);
    config_destroy(cfg);
    
    return 0;
}

