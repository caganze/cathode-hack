/**
 * @file generate_test_data.c
 * @brief Generate synthetic test data for stream detection
 * 
 * Creates CSV files with background stars and injected streams.
 * Supports RV/distance and measurement errors.
 * 
 * Usage: 
 *   ./generate_test_data -o output.csv -n 50000 --scenario basic
 *   ./generate_test_data -o output.csv --scenario rv --error 0.1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

/* Random number generator state */
static unsigned long rand_state = 12345;

static double rand_uniform(void) {
    rand_state = rand_state * 1103515245 + 12345;
    return (double)(rand_state % 1000000) / 1000000.0;
}

static double rand_gaussian(void) {
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    if (u1 < 1e-10) u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static void rand_seed(unsigned long seed) {
    rand_state = seed;
}

/* Apply fractional error to a value */
static double apply_error(double value, double error_frac) {
    if (error_frac <= 0) return value;
    return value * (1.0 + rand_gaussian() * error_frac);
}

/* Apply absolute error to a value */
static double apply_abs_error(double value, double error_abs) {
    if (error_abs <= 0) return value;
    return value + rand_gaussian() * error_abs;
}

/* Convert equatorial to galactic */
static void equatorial_to_galactic(double ra, double dec, double *l, double *b) {
    double gnp_ra = 192.85948 * DEG2RAD;
    double gnp_dec = 27.12825 * DEG2RAD;
    double l_ncp = 122.932 * DEG2RAD;
    
    double ra_rad = ra * DEG2RAD;
    double dec_rad = dec * DEG2RAD;
    
    double sin_b = sin(gnp_dec) * sin(dec_rad) + 
                   cos(gnp_dec) * cos(dec_rad) * cos(ra_rad - gnp_ra);
    *b = asin(sin_b) * RAD2DEG;
    
    double y = cos(dec_rad) * sin(ra_rad - gnp_ra);
    double x = cos(gnp_dec) * sin(dec_rad) - 
               sin(gnp_dec) * cos(dec_rad) * cos(ra_rad - gnp_ra);
    *l = fmod((l_ncp - atan2(y, x)) * RAD2DEG + 360.0, 360.0);
}

/* Convert stream coordinates to equatorial */
static void stream_to_equatorial(double pole_ra, double pole_dec, 
                                  double phi1, double phi2,
                                  double *ra, double *dec) {
    double pole_ra_rad = pole_ra * DEG2RAD;
    double pole_dec_rad = pole_dec * DEG2RAD;
    double phi1_rad = phi1 * DEG2RAD;
    double phi2_rad = phi2 * DEG2RAD;
    
    double x = cos(phi2_rad) * cos(phi1_rad);
    double y = cos(phi2_rad) * sin(phi1_rad);
    double z = sin(phi2_rad);
    
    double sin_dec_pole = sin(pole_dec_rad);
    double cos_dec_pole = cos(pole_dec_rad);
    double sin_ra_pole = sin(pole_ra_rad);
    double cos_ra_pole = cos(pole_ra_rad);
    
    double x_eq = -sin_ra_pole * y + cos_ra_pole * (cos_dec_pole * x + sin_dec_pole * z);
    double y_eq = cos_ra_pole * y + sin_ra_pole * (cos_dec_pole * x + sin_dec_pole * z);
    double z_eq = -sin_dec_pole * x + cos_dec_pole * z;
    
    if (z_eq > 1.0) z_eq = 1.0;
    if (z_eq < -1.0) z_eq = -1.0;
    
    *dec = asin(z_eq) * RAD2DEG;
    *ra = fmod(atan2(y_eq, x_eq) * RAD2DEG + 360.0, 360.0);
}

/* Extended stream parameters with RV and distance */
typedef struct {
    char name[64];
    double pole_ra, pole_dec;
    double phi1_start, phi1_end;
    double width;           /* degrees */
    double pm_ra, pm_dec;   /* mas/yr */
    double pm_dispersion;   /* mas/yr */
    double rv;              /* km/s */
    double rv_dispersion;   /* km/s */
    double distance;        /* kpc */
    double dist_dispersion; /* kpc */
    int n_stars;
} StreamParams;

/* Global error fraction */
static double g_error_frac = 0.0;
static int g_use_rv_dist = 0;

/* Background patch parameters (global for simplicity) */
/* Default: high galactic latitude region centered at (RA=180, Dec=50) */
/* This corresponds to l~270, b~+65 - well above galactic plane */
static double g_patch_center_ra = 180.0;   /* Center RA of background patch */
static double g_patch_center_dec = 50.0;   /* Center Dec of background patch */
static double g_patch_radius = 20.0;       /* Radius in degrees */

/* Generate background star within a patch */
static void generate_background_star(FILE *f, double min_b) {
    double ra, dec, l, b;
    
    /* Generate position within patch, avoiding galactic plane */
    do {
        /* Generate uniform within circular patch */
        double r = sqrt(rand_uniform()) * g_patch_radius;
        double theta = rand_uniform() * 2.0 * M_PI;
        double dra = r * cos(theta) / cos(g_patch_center_dec * DEG2RAD);
        double ddec = r * sin(theta);
        ra = g_patch_center_ra + dra;
        dec = g_patch_center_dec + ddec;
        
        /* Wrap RA */
        if (ra < 0) ra += 360.0;
        if (ra >= 360) ra -= 360.0;
        
        /* Clamp Dec */
        if (dec < -90) dec = -90;
        if (dec > 90) dec = 90;
        
        equatorial_to_galactic(ra, dec, &l, &b);
    } while (fabs(b) < min_b);
    
    /* Proper motions */
    double pmra = rand_gaussian() * 5.0;
    double pmdec = rand_gaussian() * 5.0;
    if (fabs(b) < 40) {
        pmra += rand_gaussian() * 3.0;
        pmdec += rand_gaussian() * 2.0;
    }
    
    /* RV: halo-like distribution */
    double rv = rand_gaussian() * 100.0;  /* km/s, large dispersion */
    
    /* Distance: uniform in volume out to ~20 kpc */
    double dist = pow(rand_uniform(), 1.0/3.0) * 20.0;  /* kpc */
    double parallax = 1.0 / dist;  /* mas */
    
    /* Apply errors */
    ra = apply_abs_error(ra, g_error_frac * 0.001);  /* ~1 mas position error */
    dec = apply_abs_error(dec, g_error_frac * 0.001);
    pmra = apply_error(pmra, g_error_frac);
    pmdec = apply_error(pmdec, g_error_frac);
    rv = apply_abs_error(rv, g_error_frac * 10.0);  /* ~10 km/s base error */
    parallax = apply_error(parallax, g_error_frac);
    
    if (g_use_rv_dist) {
        fprintf(f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,0\n", 
                ra, dec, pmra, pmdec, rv, parallax);
    } else {
        fprintf(f, "%.8f,%.8f,%.4f,%.4f,0\n", ra, dec, pmra, pmdec);
    }
}

/* Minimum galactic latitude for valid footprint */
#define MIN_GALACTIC_B 20.0

/* Angular distance in degrees */
static double angular_distance(double ra1, double dec1, double ra2, double dec2) {
    double ra1_rad = ra1 * DEG2RAD;
    double dec1_rad = dec1 * DEG2RAD;
    double ra2_rad = ra2 * DEG2RAD;
    double dec2_rad = dec2 * DEG2RAD;
    
    double cos_d = sin(dec1_rad) * sin(dec2_rad) + 
                   cos(dec1_rad) * cos(dec2_rad) * cos(ra1_rad - ra2_rad);
    if (cos_d > 1.0) cos_d = 1.0;
    if (cos_d < -1.0) cos_d = -1.0;
    return acos(cos_d) * RAD2DEG;
}

/* Generate stream star - returns 1 if star was written, 0 if skipped */
static int generate_stream_star(FILE *f, FILE *truth_f, StreamParams *s) {
    double phi1 = s->phi1_start + rand_uniform() * (s->phi1_end - s->phi1_start);
    double phi2 = rand_gaussian() * s->width;
    
    double ra, dec;
    stream_to_equatorial(s->pole_ra, s->pole_dec, phi1, phi2, &ra, &dec);
    
    /* Check if within valid footprint (|b| > MIN_GALACTIC_B) */
    double l, b;
    equatorial_to_galactic(ra, dec, &l, &b);
    if (fabs(b) < MIN_GALACTIC_B) {
        return 0;  /* Skip stars in galactic plane */
    }
    
    /* Check if within background patch */
    double dist = angular_distance(ra, dec, g_patch_center_ra, g_patch_center_dec);
    if (dist > g_patch_radius) {
        return 0;  /* Skip stars outside patch */
    }
    
    /* True values */
    double pmra_true = s->pm_ra + rand_gaussian() * s->pm_dispersion;
    double pmdec_true = s->pm_dec + rand_gaussian() * s->pm_dispersion;
    double rv_true = s->rv + rand_gaussian() * s->rv_dispersion;
    double dist_true = s->distance + rand_gaussian() * s->dist_dispersion;
    if (dist_true < 0.1) dist_true = 0.1;
    double parallax_true = 1.0 / dist_true;
    
    /* Apply errors */
    double ra_obs = apply_abs_error(ra, g_error_frac * 0.001);
    double dec_obs = apply_abs_error(dec, g_error_frac * 0.001);
    double pmra_obs = apply_error(pmra_true, g_error_frac);
    double pmdec_obs = apply_error(pmdec_true, g_error_frac);
    double rv_obs = apply_abs_error(rv_true, g_error_frac * 5.0);
    double parallax_obs = apply_error(parallax_true, g_error_frac);
    
    if (g_use_rv_dist) {
        fprintf(f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,1\n",
                ra_obs, dec_obs, pmra_obs, pmdec_obs, rv_obs, parallax_obs);
    } else {
        fprintf(f, "%.8f,%.8f,%.4f,%.4f,1\n", ra_obs, dec_obs, pmra_obs, pmdec_obs);
    }
    
    if (truth_f) {
        if (g_use_rv_dist) {
            fprintf(truth_f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,1,%s\n",
                    ra, dec, pmra_true, pmdec_true, rv_true, parallax_true, s->name);
        } else {
            fprintf(truth_f, "%.8f,%.8f,%.4f,%.4f,1,%s\n",
                    ra, dec, pmra_true, pmdec_true, s->name);
        }
    }
    return 1;
}

/* Generate background stars */
static int generate_background(FILE *f, int n_stars, double min_b) {
    printf("Generating %d background stars...\n", n_stars);
    
    for (int i = 0; i < n_stars; i++) {
        generate_background_star(f, min_b);
        if ((i + 1) % 10000 == 0) {
            printf("  Generated %d stars...\n", i + 1);
        }
    }
    
    return n_stars;
}

/* Generate stream stars */
static int generate_stream(FILE *f, FILE *truth_f, StreamParams *s) {
    /* Show approximate stream location */
    double mid_phi1 = (s->phi1_start + s->phi1_end) / 2.0;
    double mid_ra, mid_dec, mid_l, mid_b;
    stream_to_equatorial(s->pole_ra, s->pole_dec, mid_phi1, 0.0, &mid_ra, &mid_dec);
    equatorial_to_galactic(mid_ra, mid_dec, &mid_l, &mid_b);
    printf("Generating stream '%s': %d stars, center at (RA=%.1f, Dec=%.1f, b=%.1f)\n", 
           s->name, s->n_stars, mid_ra, mid_dec, mid_b);
    
    int generated = 0;
    int attempts = 0;
    int max_attempts = s->n_stars * 20;  /* More attempts for tricky geometries */
    
    while (generated < s->n_stars && attempts < max_attempts) {
        if (generate_stream_star(f, truth_f, s)) {
            generated++;
        }
        attempts++;
    }
    
    if (generated < s->n_stars) {
        printf("  Warning: Only generated %d/%d stars (|b|>20 constraint)\n", 
               generated, s->n_stars);
    }
    
    return generated;
}

/* Predefined streams - poles chosen so streams pass through default patch (180, 50) */
/* Stream path is great circle 90° from pole. Pole at (90, 0) gives stream through (180, ~50) */

/* Stream star counts scaled for 300 stars/sq.deg background (12x original) */

/* GD-1 like: cold, thin stream */
static StreamParams gd1_like = {
    "GD1-like", 90.0, -10.0, -20.0, 20.0, 0.3,  /* pole at RA=90, Dec=-10 */
    -12.0, -3.0, 0.5,      /* PM: mean_ra, mean_dec, dispersion */
    -150.0, 20.0,          /* RV: mean, dispersion (km/s) */
    10.0, 1.0,             /* Distance: mean, dispersion (kpc) */
    6000                   /* 500 * 12 */
};

/* Faint narrow: different orientation */
static StreamParams faint_narrow = {
    "Faint-narrow", 270.0, 10.0, -15.0, 15.0, 0.2,
    -5.0, 2.0, 0.3,
    50.0, 15.0,
    15.0, 2.0,
    2400                   /* 200 * 12 */
};

/* Bright wide: another orientation */
static StreamParams bright_wide = {
    "Bright-wide", 180.0, -40.0, -30.0, 30.0, 1.0,
    8.0, -1.0, 0.8,
    -80.0, 25.0,
    5.0, 0.5,
    9600                   /* 800 * 12 */
};

/* High PM: nearby stream with high proper motion */
static StreamParams high_pm = {
    "High-PM", 0.0, -20.0, -25.0, 25.0, 0.4,
    -18.0, 5.0, 0.4,
    200.0, 30.0,
    3.0, 0.3,
    4200                   /* 350 * 12 */
};

/* Diffuse stream - 5x wider than GD-1 (harder to detect) */
static StreamParams diffuse_5x = {
    "Diffuse-5x", 90.0, -10.0, -20.0, 20.0, 1.5,  /* Same pole as GD1, 5x width */
    -12.0, -3.0, 2.5,      /* PM dispersion: 0.5 -> 2.5 mas/yr */
    -150.0, 40.0,          /* RV dispersion also larger */
    10.0, 2.0,             /* Distance spread larger */
    9600                   /* 800 * 12 */
};

/* Very diffuse stream - 10x wider than GD-1 (very hard to detect) */
static StreamParams diffuse_10x = {
    "Diffuse-10x", 90.0, -10.0, -20.0, 20.0, 3.0,  /* Same pole, 10x width */
    -12.0, -3.0, 5.0,      /* PM dispersion: 0.5 -> 5.0 mas/yr */
    -150.0, 60.0,          /* RV dispersion even larger */
    10.0, 3.0,             /* Distance spread larger */
    14400                  /* 1200 * 12 */
};

/* ============================================================================
 * Ultra-cold streams with very low RV dispersion (< 5 km/s)
 * These are the most physically realistic streams for testing RV filtering
 * ============================================================================ */

/* Pal 5 like: very cold tidal stream from dwarf/GC disruption */
static StreamParams pal5_like = {
    "Pal5-like", 120.0, 0.0, -10.0, 10.0, 0.4,  /* Different pole */
    -2.5, -1.8, 0.3,       /* PM: low dispersion */
    -58.0, 2.0,            /* RV: mean, VERY low dispersion (2 km/s!) */
    23.0, 0.5,             /* Distance: well-constrained */
    3600                   /* 300 * 12 */
};

/* Ultra-cold stream 1: dynamically cold, recently disrupted */
static StreamParams ultra_cold_1 = {
    "UltraCold-1", 60.0, -30.0, -15.0, 15.0, 0.25,
    -8.0, 2.5, 0.2,        /* Very tight PM */
    -120.0, 3.0,           /* RV dispersion: 3 km/s */
    12.0, 0.3,             /* Tight distance */
    4800                   /* 400 * 12 */
};

/* Ultra-cold stream 2: different kinematics */
static StreamParams ultra_cold_2 = {
    "UltraCold-2", 200.0, 20.0, -20.0, 20.0, 0.35,
    5.0, -4.0, 0.25,
    80.0, 4.0,             /* RV dispersion: 4 km/s */
    8.0, 0.4,
    4200                   /* 350 * 12 */
};

/* Hot contaminated stream: good PM but high RV dispersion (should be rejected) */
static StreamParams hot_contaminant = {
    "HotContam", 90.0, -10.0, -20.0, 20.0, 0.3,  /* Same geometry as GD1 */
    -12.0, -3.0, 0.5,      /* Good PM (same as GD1) */
    0.0, 80.0,             /* HIGH RV dispersion (80 km/s!) - field contamination */
    10.0, 5.0,             /* Spread in distance */
    6000                   /* Same star count as GD1 */
};

/* Globular cluster - compact, spherical, coherent motion */
typedef struct {
    char name[64];
    double ra, dec;         /* Center position */
    double radius;          /* Angular radius (degrees) */
    double pm_ra, pm_dec;   /* Mean proper motion */
    double pm_dispersion;   /* Internal velocity dispersion */
    double rv;              /* Mean RV */
    double rv_dispersion;   /* RV dispersion */
    double distance;        /* Distance (kpc) */
    double dist_spread;     /* Distance spread */
    int n_stars;
} ClusterParams;

/* Globular clusters - positioned within default patch (RA~180, Dec~50) */
/* Star counts scaled for 300 stars/sq.deg background (12x original) */
static ClusterParams omega_cen = {
    "GC-North", 178.0, 52.0, 0.5,   /* Near patch center */
    -2.5, -1.5, 0.8,
    -50.0, 15.0,
    12.0, 0.5,
    12000                  /* 1000 * 12 */
};

static ClusterParams m4_like = {
    "GC-South", 182.0, 48.0, 0.3,   /* Also near patch center */
    1.0, -2.0, 0.5,
    80.0, 10.0,
    8.0, 0.3,
    6000                   /* 500 * 12 */
};

/* Generate single cluster star - returns 1 if written, 0 if skipped */
static int generate_cluster_star(FILE *f, FILE *truth_f, ClusterParams *c) {
    /* Position: Gaussian distribution around center */
    double r = fabs(rand_gaussian()) * c->radius / 2.0;
    double theta = rand_uniform() * 2.0 * M_PI;
    double ra = c->ra + r * cos(theta) / cos(c->dec * DEG2RAD);
    double dec = c->dec + r * sin(theta);
    
    /* Check if within valid footprint */
    double l, b;
    equatorial_to_galactic(ra, dec, &l, &b);
    if (fabs(b) < MIN_GALACTIC_B) {
        return 0;
    }
    
    /* Check if within background patch */
    double dist = angular_distance(ra, dec, g_patch_center_ra, g_patch_center_dec);
    if (dist > g_patch_radius) {
        return 0;
    }
    
    /* True values with internal dispersion */
    double pmra_true = c->pm_ra + rand_gaussian() * c->pm_dispersion;
    double pmdec_true = c->pm_dec + rand_gaussian() * c->pm_dispersion;
    double rv_true = c->rv + rand_gaussian() * c->rv_dispersion;
    double dist_true = c->distance + rand_gaussian() * c->dist_spread;
    if (dist_true < 0.1) dist_true = 0.1;
    double parallax_true = 1.0 / dist_true;
    
    /* Apply measurement errors */
    double ra_obs = apply_abs_error(ra, g_error_frac * 0.001);
    double dec_obs = apply_abs_error(dec, g_error_frac * 0.001);
    double pmra_obs = apply_error(pmra_true, g_error_frac);
    double pmdec_obs = apply_error(pmdec_true, g_error_frac);
    double rv_obs = apply_abs_error(rv_true, g_error_frac * 5.0);
    double parallax_obs = apply_error(parallax_true, g_error_frac);
    
    if (g_use_rv_dist) {
        fprintf(f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,1\n",
                ra_obs, dec_obs, pmra_obs, pmdec_obs, rv_obs, parallax_obs);
    } else {
        fprintf(f, "%.8f,%.8f,%.4f,%.4f,1\n", ra_obs, dec_obs, pmra_obs, pmdec_obs);
    }
    
    if (truth_f) {
        if (g_use_rv_dist) {
            fprintf(truth_f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,1,%s\n",
                    ra, dec, pmra_true, pmdec_true, rv_true, parallax_true, c->name);
        } else {
            fprintf(truth_f, "%.8f,%.8f,%.4f,%.4f,1,%s\n",
                    ra, dec, pmra_true, pmdec_true, c->name);
        }
    }
    return 1;
}

/* Generate globular cluster stars */
static int generate_cluster(FILE *f, FILE *truth_f, ClusterParams *c) {
    /* First check if cluster center is in valid footprint */
    double l, b;
    equatorial_to_galactic(c->ra, c->dec, &l, &b);
    printf("Generating cluster '%s' at (RA=%.1f, Dec=%.1f, b=%.1f) with %d stars...\n", 
           c->name, c->ra, c->dec, b, c->n_stars);
    
    if (fabs(b) < MIN_GALACTIC_B) {
        printf("  Warning: Cluster center at |b|=%.1f is in galactic plane, skipping!\n", fabs(b));
        return 0;
    }
    
    int generated = 0;
    int attempts = 0;
    int max_attempts = c->n_stars * 10;
    
    while (generated < c->n_stars && attempts < max_attempts) {
        if (generate_cluster_star(f, truth_f, c)) {
            generated++;
        }
        attempts++;
    }
    
    if (generated < c->n_stars) {
        printf("  Warning: Only generated %d/%d stars (footprint constraint)\n",
               generated, c->n_stars);
    }
    
    return generated;
}

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  -o, --output FILE    Output CSV file (default: test_data.csv)\n");
    printf("  -d, --density D      Background density in stars/sq.deg (default: 300)\n");
    printf("  -n, --n-background N Override: exact number of background stars\n");
    printf("  -r, --patch-radius R Background patch radius in deg (default: 20)\n");
    printf("  --patch-ra RA        Background patch center RA (default: 180)\n");
    printf("  --patch-dec DEC      Background patch center Dec (default: 45)\n");
    printf("  -s, --seed SEED      Random seed (default: 42)\n");
    printf("  --scenario NAME      Test scenario:\n");
    printf("                         basic    - Single GD1-like stream\n");
    printf("                         gd1      - GD1-like stream (1000 stars)\n");
    printf("                         diffuse5 - 5x more diffuse than GD1\n");
    printf("                         diffuse10- 10x more diffuse than GD1\n");
    printf("                         gc       - Globular cluster (M53-like)\n");
    printf("                         multi    - Multiple streams\n");
    printf("                         gc_multi - Streams + globular clusters\n");
    printf("                         rv       - Include RV and distance\n");
    printf("                         rv_multi - RV/distance with multiple streams\n");
    printf("                         rv_gc    - RV/distance with globular cluster\n");
    printf("                         rv_diffuse - RV/distance with diffuse streams\n");
    printf("                         rv_cold  - Ultra-cold streams (disp < 5 km/s)\n");
    printf("                         rv_mixed - Cold streams + hot field contaminants\n");
    printf("  --error FRAC         Measurement error fraction (0.0, 0.1, 0.5)\n");
    printf("  -h, --help           Show this help\n");
    printf("\nBackground density:\n");
    printf("  Default: 300 stars/sq.deg with 20 deg radius = ~377,000 stars\n");
    printf("  Use -r to change patch size: -r 10 gives ~94,000 stars\n");
    printf("\nExamples:\n");
    printf("  %s --scenario basic                  # 300/sq.deg, no errors\n", prog);
    printf("  %s --scenario gd1 -r 10              # Smaller patch (94k stars)\n", prog);
    printf("  %s --scenario rv --error 0.1         # RV+dist, 10%% errors\n", prog);
}

int main(int argc, char *argv[]) {
    char output[256] = "test_data.csv";
    int n_background = 50000;
    unsigned long seed = 42;
    char scenario[32] = "basic";
    double error_frac = 0.0;
    
    static struct option long_options[] = {
        {"output", required_argument, 0, 'o'},
        {"n-background", required_argument, 0, 'n'},
        {"density", required_argument, 0, 'd'},
        {"patch-radius", required_argument, 0, 'r'},
        {"patch-ra", required_argument, 0, 'R'},
        {"patch-dec", required_argument, 0, 'D'},
        {"seed", required_argument, 0, 's'},
        {"scenario", required_argument, 0, 'S'},
        {"error", required_argument, 0, 'e'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    double density = 300.0;  /* Default: 300 stars/sq.deg */
    int use_density = 1;     /* Default: use density mode */
    
    int opt;
    while ((opt = getopt_long(argc, argv, "o:n:d:r:s:e:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o': strncpy(output, optarg, sizeof(output)-1); break;
            case 'n': n_background = atoi(optarg); use_density = 0; break;
            case 'd': density = atof(optarg); use_density = 1; break;
            case 'r': g_patch_radius = atof(optarg); break;
            case 'R': g_patch_center_ra = atof(optarg); break;
            case 'D': g_patch_center_dec = atof(optarg); break;
            case 's': seed = atol(optarg); break;
            case 'S': strncpy(scenario, optarg, sizeof(scenario)-1); break;
            case 'e': error_frac = atof(optarg); break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }
    
    /* Set global error */
    g_error_frac = error_frac;
    
    /* Calculate n_background from density if using density mode */
    double patch_area = M_PI * g_patch_radius * g_patch_radius;  /* sq.deg */
    if (use_density) {
        n_background = (int)(density * patch_area);
    }
    
    /* Check if scenario uses RV/distance */
    g_use_rv_dist = (strstr(scenario, "rv") != NULL);
    
    printf("Synthetic Data Generator\n");
    printf("========================\n");
    printf("Output: %s\n", output);
    printf("Patch: (RA=%.1f, Dec=%.1f), radius=%.1f deg, area=%.1f sq.deg\n",
           g_patch_center_ra, g_patch_center_dec, g_patch_radius, patch_area);
    printf("Density: %.0f stars/sq.deg\n", (double)n_background / patch_area);
    printf("Background: %d stars\n", n_background);
    printf("Scenario: %s\n", scenario);
    printf("Error fraction: %.0f%%\n", error_frac * 100);
    printf("Include RV/dist: %s\n", g_use_rv_dist ? "yes" : "no");
    printf("Seed: %lu\n\n", seed);
    
    rand_seed(seed);
    
    /* Open output file */
    FILE *f = fopen(output, "w");
    if (!f) {
        fprintf(stderr, "Error: Cannot open %s for writing\n", output);
        return 1;
    }
    
    /* Open truth file */
    char truth_file[256];
    snprintf(truth_file, sizeof(truth_file), "%s", output);
    char *ext = strrchr(truth_file, '.');
    if (ext) strcpy(ext, "_truth.csv");
    else strcat(truth_file, "_truth");
    
    FILE *truth_f = fopen(truth_file, "w");
    
    /* Write headers */
    if (g_use_rv_dist) {
        fprintf(f, "ra,dec,pmra,pmdec,rv,parallax,is_stream\n");
        if (truth_f) fprintf(truth_f, "ra,dec,pmra,pmdec,rv,parallax,is_stream,stream_name\n");
    } else {
        fprintf(f, "ra,dec,pmra,pmdec,is_stream\n");
        if (truth_f) fprintf(truth_f, "ra,dec,pmra,pmdec,is_stream,stream_name\n");
    }
    
    /* Generate background */
    int total = generate_background(f, n_background, 20.0);
    
    /* Generate streams and clusters based on scenario */
    int n_streams = 0;
    StreamParams *streams[10];
    int n_clusters = 0;
    ClusterParams *clusters[10];
    
    if (strcmp(scenario, "basic") == 0) {
        streams[n_streams++] = &gd1_like;
    } else if (strcmp(scenario, "gd1") == 0) {
        gd1_like.n_stars = 1000;
        streams[n_streams++] = &gd1_like;
    } else if (strcmp(scenario, "diffuse5") == 0) {
        /* 5x more diffuse than GD-1 */
        streams[n_streams++] = &diffuse_5x;
    } else if (strcmp(scenario, "diffuse10") == 0) {
        /* 10x more diffuse than GD-1 */
        streams[n_streams++] = &diffuse_10x;
    } else if (strcmp(scenario, "diffuse_compare") == 0) {
        /* Compare: GD1 + 5x + 10x diffuse side by side */
        streams[n_streams++] = &gd1_like;
        streams[n_streams++] = &diffuse_5x;
        streams[n_streams++] = &diffuse_10x;
    } else if (strcmp(scenario, "gc") == 0) {
        clusters[n_clusters++] = &omega_cen;
    } else if (strcmp(scenario, "multi") == 0) {
        streams[n_streams++] = &gd1_like;
        streams[n_streams++] = &faint_narrow;
        streams[n_streams++] = &bright_wide;
        streams[n_streams++] = &high_pm;
    } else if (strcmp(scenario, "gc_multi") == 0) {
        streams[n_streams++] = &gd1_like;
        streams[n_streams++] = &faint_narrow;
        clusters[n_clusters++] = &omega_cen;
        clusters[n_clusters++] = &m4_like;
    } else if (strcmp(scenario, "rv") == 0) {
        streams[n_streams++] = &gd1_like;
    } else if (strcmp(scenario, "rv_multi") == 0) {
        streams[n_streams++] = &gd1_like;
        streams[n_streams++] = &faint_narrow;
        streams[n_streams++] = &bright_wide;
        streams[n_streams++] = &high_pm;
    } else if (strcmp(scenario, "rv_gc") == 0) {
        clusters[n_clusters++] = &omega_cen;
    } else if (strcmp(scenario, "rv_diffuse") == 0) {
        /* Diffuse streams with RV/distance */
        streams[n_streams++] = &diffuse_5x;
        streams[n_streams++] = &diffuse_10x;
    } else if (strcmp(scenario, "rv_cold") == 0) {
        /* Ultra-cold streams with very low RV dispersion (< 5 km/s) */
        /* These should be easily detectable with RV filtering */
        streams[n_streams++] = &pal5_like;
        streams[n_streams++] = &ultra_cold_1;
        streams[n_streams++] = &ultra_cold_2;
    } else if (strcmp(scenario, "rv_mixed") == 0) {
        /* Mix of cold streams and hot contaminants */
        /* Tests that RV dispersion filtering rejects hot structures */
        streams[n_streams++] = &pal5_like;      /* Cold: should be detected */
        streams[n_streams++] = &ultra_cold_1;   /* Cold: should be detected */
        streams[n_streams++] = &hot_contaminant;/* Hot: should be rejected */
    } else if (strcmp(scenario, "rv_pal5") == 0) {
        /* Single Pal 5-like stream for focused testing */
        streams[n_streams++] = &pal5_like;
    } else if (strcmp(scenario, "rv_hot") == 0) {
        /* Only hot contaminant - should be fully rejected */
        streams[n_streams++] = &hot_contaminant;
    } else {
        fprintf(stderr, "Unknown scenario: %s\n", scenario);
        fclose(f);
        if (truth_f) fclose(truth_f);
        return 1;
    }
    
    int stream_stars = 0;
    for (int i = 0; i < n_streams; i++) {
        stream_stars += generate_stream(f, truth_f, streams[i]);
    }
    
    int cluster_stars = 0;
    for (int i = 0; i < n_clusters; i++) {
        cluster_stars += generate_cluster(f, truth_f, clusters[i]);
    }
    
    total += stream_stars + cluster_stars;
    
    fclose(f);
    if (truth_f) fclose(truth_f);
    
    printf("\nGenerated %d total stars\n", total);
    printf("  Background: %d\n", n_background);
    printf("  Stream: %d\n", stream_stars);
    printf("  Cluster: %d\n", cluster_stars);
    printf("\nFiles written:\n");
    printf("  %s\n", output);
    printf("  %s\n", truth_file);
    
    return 0;
}
