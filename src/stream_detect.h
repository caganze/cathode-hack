/**
 * @file stream_detect.h
 * @brief Via Machinae-inspired stellar stream detection using anomaly detection
 * 
 * Implementation based on arXiv:2509.08064v1
 * "Via Machinae 3.0: A search for stellar streams in Gaia with the CATHODE algorithm"
 * 
 * Adapted for DESI DR1 Milky Way Survey catalog data.
 * 
 * NOTE: This implementation uses only positions (RA, Dec) and proper motions 
 * (pmra, pmdec) for stream detection. Isochrone-based cuts (magnitude, color, 
 * f_dim) have been removed for simplicity.
 */

#ifndef STREAM_DETECT_H
#define STREAM_DETECT_H

#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Constants and Configuration
 * ============================================================================ */

#define MAX_STARS_PER_PATCH    100000    /* Maximum stars in a single patch */
#define MAX_PATCHES            5000      /* Maximum number of sky patches (more small patches) */
#define MAX_ROIS_PER_PATCH     1000      /* Maximum ROIs per patch */
#define MAX_STARS_PER_ROI      200       /* Stars selected per ROI */
#define MAX_LINE_STARS         100       /* Stars for line fitting */
#define MAX_PROTOCLUSTERS      1000      /* Maximum protoclusters */
#define MAX_PROTOSTREAMS       500       /* Maximum protostreams */
#define MAX_STREAM_CANDIDATES  5000      /* Maximum final stream candidates */

/* Algorithm parameters */
/* Sky patches: 10 sq.deg -> radius = sqrt(10/pi) ~ 1.78 deg */
#define PATCH_RADIUS_DEG       1.78      /* Patch radius for ~10 sq.deg area */
#define PATCH_AREA_SQDEG       10.0      /* Target patch area in sq.deg */

/* Proper motion signal regions: 0.5 mas/yr width */
#define SR_WIDTH_MAS_YR        0.5       /* Signal region width in mas/yr */
#define SR_STEP_MAS_YR         0.25      /* Signal region step in mas/yr */
#define ROI_WIDTH_MAS_YR       0.5       /* ROI width in mas/yr */

/* Other parameters */
#define FIDUCIAL_DIST_DEG      1.78      /* Same as patch radius */
#define MIN_STARS_IN_SR        10        /* Minimum stars in signal region */
#define MAX_STARS_IN_SR        100000    /* Maximum stars in signal region */
#define TOP_N_ANOMALOUS        30        /* Number of most anomalous stars to select */
#define PROTOCLUSTER_SIG_CUT   5.0       /* Protocluster significance cutoff */
#define EDGE_RADIUS_DEG        1.5       /* Edge protostream radius threshold */
#define DUPLICATE_OVERLAP      0.4       /* Overlap threshold for duplicates */
/* Isochrone cuts removed - algorithm uses only positions and proper motions */

/* Hough transform parameters */
#define HOUGH_THETA_BINS       180       /* Number of theta bins */
#define HOUGH_RHO_BINS         200       /* Number of rho bins */
#define HOUGH_THETA_MIN        -M_PI     /* Minimum theta */
#define HOUGH_THETA_MAX        M_PI      /* Maximum theta */

/* KDE bandwidth (Silverman's rule will be used to compute optimal) */
#define KDE_DEFAULT_BANDWIDTH  0.5

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/**
 * @brief Single star data structure
 */
typedef struct {
    /* Original equatorial coordinates */
    double ra;              /* Right ascension (degrees) */
    double dec;             /* Declination (degrees) */
    
    /* Proper motions */
    double pmra;            /* Proper motion in RA (mas/yr) */
    double pmdec;           /* Proper motion in Dec (mas/yr) */
    
    /* Photometry */
    double g_mag;           /* G-band magnitude */
    double bp_rp;           /* BP-RP color (or similar) */
    
    /* Additional from DESI */
    double rv;              /* Radial velocity (km/s) */
    double rv_err;          /* Radial velocity error */
    double parallax;        /* Parallax (mas) */
    double parallax_err;    /* Parallax error (mas) */
    double distance;        /* Distance (kpc, computed from parallax or rvsdistnn) */
    double teff;            /* Effective temperature (K) */
    double teff_err;        /* Teff uncertainty */
    double logg;            /* Surface gravity */
    double logg_err;        /* Logg uncertainty */
    double feh;             /* Metallicity [Fe/H] */
    double feh_err;         /* [Fe/H] uncertainty */
    double vsini;           /* Rotational velocity (km/s) */
    
    /* Patch-local coordinates (computed) */
    double phi;             /* Local longitude (degrees) */
    double lambda;          /* Local latitude (degrees) */
    double mu_phi;          /* Proper motion in phi (mas/yr) */
    double mu_lambda;       /* Proper motion in lambda (mas/yr) */
    
    /* Computed quantities */
    double anomaly_score;   /* R(x) anomaly score */
    double dist_from_center;/* Angular distance from patch center */
    
    /* Identifiers */
    int64_t target_id;      /* DESI target ID */
    uint32_t star_idx;      /* Index in original array */
    
    /* Flags */
    bool is_selected;       /* Selected for further analysis */
    bool passes_fiducial;   /* Passes fiducial cuts */
} Star;

/**
 * @brief Sky patch containing stars
 */
typedef struct {
    /* Patch center */
    double center_ra;       /* Center RA (degrees) */
    double center_dec;      /* Center Dec (degrees) */
    double center_l;        /* Center galactic longitude */
    double center_b;        /* Center galactic latitude */
    
    /* Patch geometry */
    double radius;          /* Patch radius (degrees) */
    
    /* Stars in this patch */
    Star *stars;            /* Array of stars */
    uint32_t n_stars;       /* Number of stars */
    uint32_t capacity;      /* Allocated capacity */
    
    /* Patch identification */
    int patch_id;           /* Unique patch ID */
    bool is_valid;          /* Patch passes quality cuts */
} Patch;

/**
 * @brief Signal Region definition
 */
typedef struct {
    double mu_min;          /* Minimum proper motion */
    double mu_max;          /* Maximum proper motion */
    bool use_mu_lambda;     /* true=mu_lambda, false=mu_phi */
    
    Star **sr_stars;        /* Pointers to stars in SR */
    uint32_t n_sr_stars;    /* Number of stars in SR */
    
    Star **sb_stars;        /* Pointers to stars in SB */
    uint32_t n_sb_stars;    /* Number of stars in SB */
} SignalRegion;

/**
 * @brief Region of Interest within a Signal Region
 */
typedef struct {
    double mu_perp_min;     /* Minimum orthogonal proper motion */
    double mu_perp_max;     /* Maximum orthogonal proper motion */
    
    Star **roi_stars;       /* Pointers to stars in ROI */
    uint32_t n_roi_stars;   /* Number of stars in ROI */
    
    /* Selected anomalous stars */
    Star **selected_stars;  /* Top anomalous stars */
    uint32_t n_selected;    /* Number of selected stars */
    
    /* Line fitting results */
    double theta_line;      /* Best-fit line angle */
    double rho_line;        /* Best-fit line distance from center */
    double sigma_line;      /* Line significance */
    
    /* Identification */
    int roi_id;             /* ROI ID within SR */
    int sr_id;              /* Parent SR ID */
    int patch_id;           /* Parent patch ID */
} ROI;

/**
 * @brief Hough accumulator for line detection
 */
typedef struct {
    int *accumulator;       /* 2D accumulator array */
    int n_theta;            /* Number of theta bins */
    int n_rho;              /* Number of rho bins */
    double theta_min;       /* Minimum theta */
    double theta_max;       /* Maximum theta */
    double rho_max;         /* Maximum rho (patch radius) */
    double d_theta;         /* Theta bin width */
    double d_rho;           /* Rho bin width */
} HoughAccumulator;

/**
 * @brief Protocluster: grouped ROIs with significant line
 */
typedef struct {
    ROI **rois;             /* Array of ROI pointers */
    uint32_t n_rois;        /* Number of ROIs */
    
    Star **line_stars;      /* Stars on the line */
    uint32_t n_line_stars;  /* Number of line stars */
    
    /* Combined line parameters */
    double theta_line;      /* Combined line angle */
    double rho_line;        /* Combined line distance */
    double significance;    /* Combined significance */
    
    /* Proper motion statistics */
    double mean_mu_phi;     /* Mean mu_phi */
    double mean_mu_lambda;  /* Mean mu_lambda */
    double std_mu_phi;      /* Std of mu_phi */
    double std_mu_lambda;   /* Std of mu_lambda */
    
    /* Identification */
    int cluster_id;         /* Cluster ID */
    int patch_id;           /* Parent patch ID */
    bool is_edge;           /* Is this an edge cluster */
} Protocluster;

/**
 * @brief Protostream: deduplicated protoclusters in a patch
 */
typedef struct {
    Protocluster **clusters;/* Array of protocluster pointers */
    uint32_t n_clusters;    /* Number of protoclusters */
    
    double significance;    /* Highest significance among clusters */
    
    /* Identification */
    int protostream_id;     /* Protostream ID */
    int patch_id;           /* Parent patch ID */
    bool is_edge;           /* Is edge protostream */
    bool passes_cuts;       /* Passes all cuts */
} Protostream;

/**
 * @brief Stream candidate: merged protostreams across patches
 */
typedef struct {
    Protostream **protostreams; /* Array of protostream pointers */
    uint32_t n_protostreams;    /* Number of protostreams */
    
    Star **all_stars;       /* All stars in candidate */
    uint32_t n_stars;       /* Total number of stars */
    
    double significance;    /* Quadrature sum of significances */
    
    /* Position summary */
    double mean_ra;         /* Mean RA */
    double mean_dec;        /* Mean Dec */
    double mean_l;          /* Mean galactic l */
    double mean_b;          /* Mean galactic b */
    
    /* Proper motion summary */
    double mean_pmra;       /* Mean proper motion RA */
    double mean_pmdec;      /* Mean proper motion Dec */
    
    /* Identification */
    int candidate_id;       /* Candidate ID */
    bool is_known;          /* Matches known stream */
    char name[64];          /* Stream name if known/assigned */
} StreamCandidate;

/**
 * @brief Kernel Density Estimator
 * 
 * Supports per-dimension bandwidths to properly smooth over
 * survey artifacts (larger bandwidth for sky position).
 */
typedef struct {
    double bandwidth;       /* Scalar bandwidth (backward compat) */
    double *bandwidths;     /* Per-dimension bandwidths (if non-NULL) */
    double *data;           /* Data points (flattened) */
    uint32_t n_points;      /* Number of data points */
    uint32_t n_dims;        /* Number of dimensions */
} KDE;

/**
 * @brief Configuration for the stream detection algorithm
 */
typedef struct {
    /* Patch parameters */
    double patch_radius;
    double min_galactic_b;  /* Exclude galactic disk |b| < this */
    
    /* Signal region parameters */
    double sr_width;
    double sr_step;
    int min_stars_sr;
    int max_stars_sr;
    
    /* ROI parameters */
    double roi_width;
    int top_n_anomalous;
    
    /* Fiducial cuts - position-based only (no isochrone) */
    double dist_limit;
    double z_min_kpc;       /* Minimum Galactocentric altitude */
    
    /* Significance cuts */
    double protocluster_sig_cut;
    double edge_radius;
    
    /* Clustering parameters */
    double duplicate_overlap;
    double line_overlap_threshold;
    double theta_merge_threshold;
    double rho_merge_threshold;
    
    /* Sky region selection (default: full sky) */
    double ra_min;          /* Minimum RA (degrees), default 0 */
    double ra_max;          /* Maximum RA (degrees), default 360 */
    double dec_min;         /* Minimum Dec (degrees), default -90 */
    double dec_max;         /* Maximum Dec (degrees), default 90 */
    
    /* Distance cuts (kpc) */
    double dist_min_kpc;    /* Minimum distance (kpc), default 0 */
    double dist_max_kpc;    /* Maximum distance (kpc), default 1000 */
    
    /* KNN dimension options */
    bool use_distance_knn;  /* Include distance in KNN density (default: false) */
    bool use_rv_knn;        /* Include RV in KNN density and clustering (default: false) */
    
    /* Density estimation smoothing scales */
    double min_sky_bandwidth;  /* Minimum bandwidth for sky position (deg, default: 2.0) */
                               /* Should be >= survey tile spacing to smooth over artifacts */
    
    /* Output */
    char output_dir[256];
    bool verbose;
} Config;

/* ============================================================================
 * Function Declarations - Core Algorithm
 * ============================================================================ */

/* Initialization and cleanup */
Config *config_create_default(void);
void config_destroy(Config *cfg);

Patch *patch_create(double center_ra, double center_dec, double radius);
void patch_destroy(Patch *patch);
int patch_add_star(Patch *patch, const Star *star);

/* Coordinate transformations */
void equatorial_to_galactic(double ra, double dec, double *l, double *b);
void galactic_to_equatorial(double l, double b, double *ra, double *dec);
void transform_to_patch_coords(Star *star, double center_ra, double center_dec);
void transform_proper_motions(Star *star, double center_ra, double center_dec);

/* Sky division */
int generate_patches(Patch **patches, int *n_patches, const Config *cfg);
int assign_stars_to_patches(Star *stars, uint32_t n_stars, 
                           Patch **patches, int n_patches);

/* Signal regions and ROIs */
SignalRegion *sr_create(Patch *patch, double mu_min, double mu_max, bool use_mu_lambda);
void sr_destroy(SignalRegion *sr);
ROI *roi_create(SignalRegion *sr, double mu_perp_min, double mu_perp_max);
void roi_destroy(ROI *roi);

/* ============================================================================
 * Function Declarations - Density Estimation
 * ============================================================================ */

KDE *kde_create(uint32_t n_dims);
void kde_destroy(KDE *kde);
int kde_fit(KDE *kde, double **data, uint32_t n_points);
double kde_evaluate(const KDE *kde, const double *point);
double kde_bandwidth_silverman(const double *data, uint32_t n_points);

/* Anomaly score computation */
int compute_anomaly_scores_kde(Patch *patch, SignalRegion *sr, const Config *cfg);
int compute_anomaly_scores_histogram(Patch *patch, SignalRegion *sr, const Config *cfg);
int compute_anomaly_scores_knn(Patch *patch, SignalRegion *sr, int k, const Config *cfg);

/* ============================================================================
 * Function Declarations - Hough Transform
 * ============================================================================ */

HoughAccumulator *hough_create(int n_theta, int n_rho, double rho_max);
void hough_destroy(HoughAccumulator *ha);
void hough_clear(HoughAccumulator *ha);
void hough_add_point(HoughAccumulator *ha, double x, double y);
int hough_find_peak(const HoughAccumulator *ha, double *theta, double *rho, int *votes);
double hough_line_significance(const HoughAccumulator *ha, int votes, int n_points);

int find_line_in_roi(ROI *roi, const Config *cfg);

/* ============================================================================
 * Function Declarations - Clustering
 * ============================================================================ */

Protocluster *protocluster_create(void);
void protocluster_destroy(Protocluster *pc);
int protocluster_add_roi(Protocluster *pc, ROI *roi);
double protocluster_compute_significance(Protocluster *pc);

Protostream *protostream_create(void);
void protostream_destroy(Protostream *ps);
int protostream_add_cluster(Protostream *ps, Protocluster *pc);

StreamCandidate *stream_candidate_create(void);
void stream_candidate_destroy(StreamCandidate *sc);
int stream_candidate_add_protostream(StreamCandidate *sc, Protostream *ps);
double stream_candidate_compute_significance(StreamCandidate *sc);

/* Clustering algorithms */
int cluster_rois_to_protoclusters(ROI **rois, int n_rois, 
                                  Protocluster **clusters, int *n_clusters,
                                  const Config *cfg);
int deduplicate_protoclusters(Protocluster **clusters, int n_clusters,
                             Protostream **protostreams, int *n_protostreams);
int merge_protostreams_across_patches(Protostream **protostreams, int n_protostreams,
                                      StreamCandidate **candidates, int *n_candidates,
                                      const Config *cfg);

/* ============================================================================
 * Function Declarations - I/O
 * ============================================================================ */

/* FITS file reading (requires cfitsio) */
int read_desi_mws_fits(const char *filename, Star **stars, uint32_t *n_stars);

/* CSV file reading (simple alternative) */
int read_stars_csv(const char *filename, Star **stars, uint32_t *n_stars);

/* Output */
int write_stream_candidates_csv(const char *filename, 
                                StreamCandidate **candidates, int n_candidates);
int write_stream_stars_csv(const char *filename, StreamCandidate *candidate);
int write_patch_summary(const char *filename, Patch **patches, int n_patches);

/* ============================================================================
 * Function Declarations - Utility
 * ============================================================================ */

/* Angular calculations */
double angular_separation(double ra1, double dec1, double ra2, double dec2);
double great_circle_distance(double ra1, double dec1, double ra2, double dec2);

/* Statistical functions */
double mean(const double *data, uint32_t n);
double variance(const double *data, uint32_t n);
double std_dev(const double *data, uint32_t n);
double median(double *data, uint32_t n);
double percentile(double *data, uint32_t n, double p);

/* Sorting */
int compare_stars_by_anomaly_score(const void *a, const void *b);
int compare_stars_by_magnitude(const void *a, const void *b);

/* Memory management */
void *safe_malloc(size_t size);
void *safe_calloc(size_t nmemb, size_t size);
void *safe_realloc(void *ptr, size_t size);

/* ============================================================================
 * Main Pipeline
 * ============================================================================ */

/**
 * @brief Run the full stream detection pipeline
 * 
 * @param stars Input star array
 * @param n_stars Number of stars
 * @param cfg Configuration
 * @param candidates Output stream candidates
 * @param n_candidates Number of candidates found
 * @return 0 on success, error code otherwise
 */
int run_stream_detection(Star *stars, uint32_t n_stars, const Config *cfg,
                         StreamCandidate ***candidates, int *n_candidates);

#endif /* STREAM_DETECT_H */

