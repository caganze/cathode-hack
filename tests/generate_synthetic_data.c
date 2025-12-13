/**
 * @file generate_synthetic_data.c
 * @brief Generate synthetic stellar data with injected streams for testing
 * 
 * Creates realistic mock catalogs with:
 * - Smooth halo/disk background stars
 * - Injected stellar streams with known parameters
 * - Various stream configurations for testing
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

/* Random number generators */
static unsigned long seed = 12345;

double rand_uniform(void) {
    seed = seed * 1103515245 + 12345;
    return (double)(seed % 1000000) / 1000000.0;
}

double rand_gaussian(void) {
    /* Box-Muller transform */
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    if (u1 < 1e-10) u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ============================================================================
 * Stream Definition
 * ============================================================================ */

typedef struct {
    char name[64];
    
    /* Stream geometry - great circle on sky */
    double pole_ra;         /* RA of stream pole (degrees) */
    double pole_dec;        /* Dec of stream pole (degrees) */
    double phi1_start;      /* Start angle along stream (degrees) */
    double phi1_end;        /* End angle along stream (degrees) */
    double width;           /* Angular width (degrees) */
    
    /* Proper motion */
    double pm_parallel;     /* PM along stream (mas/yr) */
    double pm_perp;         /* PM perpendicular to stream (mas/yr) */
    double pm_dispersion;   /* PM dispersion (mas/yr) */
    
    /* Photometry */
    double mean_mag;        /* Mean G magnitude */
    double mag_spread;      /* Magnitude spread */
    double mean_color;      /* Mean BP-RP color */
    double color_spread;    /* Color spread */
    
    /* Number of stars */
    int n_stars;
} StreamParams;

/* ============================================================================
 * Background Generation
 * ============================================================================ */

/**
 * @brief Generate smooth background stars
 * 
 * Creates a population mimicking halo + thick disk stars
 */
void generate_background(FILE *fp, int n_stars, double min_b) {
    printf("Generating %d background stars...\n", n_stars);
    
    for (int i = 0; i < n_stars; i++) {
        /* Position: uniform on sphere, avoiding Galactic plane */
        double ra, dec, l, b;
        do {
            ra = rand_uniform() * 360.0;
            /* Dec: uniform in sin(dec) */
            dec = asin(2.0 * rand_uniform() - 1.0) * RAD2DEG;
            
            /* Convert to Galactic */
            double gnp_ra = 192.85948 * DEG2RAD;
            double gnp_dec = 27.12825 * DEG2RAD;
            double l_ncp = 122.932 * DEG2RAD;
            
            double ra_rad = ra * DEG2RAD;
            double dec_rad = dec * DEG2RAD;
            
            double sin_b = sin(gnp_dec) * sin(dec_rad) + 
                          cos(gnp_dec) * cos(dec_rad) * cos(ra_rad - gnp_ra);
            b = asin(sin_b) * RAD2DEG;
            
            double y = cos(dec_rad) * sin(ra_rad - gnp_ra);
            double x = cos(gnp_dec) * sin(dec_rad) - 
                      sin(gnp_dec) * cos(dec_rad) * cos(ra_rad - gnp_ra);
            l = fmod((l_ncp - atan2(y, x)) * RAD2DEG + 360.0, 360.0);
            
        } while (fabs(b) < min_b);
        
        /* Proper motion: smooth distribution centered near 0 */
        /* Halo-like: broader distribution */
        double pmra = rand_gaussian() * 5.0;
        double pmdec = rand_gaussian() * 5.0;
        
        /* Add disk contribution for |b| < 40 */
        if (fabs(b) < 40) {
            pmra += rand_gaussian() * 3.0 + 2.0 * (rand_uniform() - 0.5);
            pmdec += rand_gaussian() * 2.0;
        }
        
        /* Magnitude: exponential distribution (more faint stars) */
        double g_mag = 15.0 + 4.0 * (-log(rand_uniform() + 0.01));
        if (g_mag > 20.5) g_mag = 20.5;
        if (g_mag < 14.0) g_mag = 14.0;
        
        /* Color: roughly solar-type main sequence */
        double bp_rp = 0.6 + rand_gaussian() * 0.3;
        if (bp_rp < 0.0) bp_rp = 0.0;
        if (bp_rp > 2.0) bp_rp = 2.0;
        
        fprintf(fp, "%.8f,%.8f,%.4f,%.4f,%.3f,%.3f,0.0,5500,4.5,-0.5,0,%d\n",
                ra, dec, pmra, pmdec, g_mag, bp_rp, i);
    }
}

/* ============================================================================
 * Stream Generation
 * ============================================================================ */

/**
 * @brief Convert stream coordinates to equatorial
 * 
 * Given a pole and angle along the stream, compute RA/Dec
 */
void stream_to_equatorial(double pole_ra, double pole_dec, 
                         double phi1, double phi2,
                         double *ra, double *dec) {
    /* phi1: angle along stream, phi2: perpendicular offset */
    double pole_ra_rad = pole_ra * DEG2RAD;
    double pole_dec_rad = pole_dec * DEG2RAD;
    double phi1_rad = phi1 * DEG2RAD;
    double phi2_rad = phi2 * DEG2RAD;
    
    /* Construct rotation matrix from stream frame to equatorial */
    /* Stream frame: z along pole, x toward phi1=0 */
    
    /* Position in stream frame */
    double x = cos(phi2_rad) * cos(phi1_rad);
    double y = cos(phi2_rad) * sin(phi1_rad);
    double z = sin(phi2_rad);
    
    /* Rotate to equatorial frame */
    double sin_dec_pole = sin(pole_dec_rad);
    double cos_dec_pole = cos(pole_dec_rad);
    double sin_ra_pole = sin(pole_ra_rad);
    double cos_ra_pole = cos(pole_ra_rad);
    
    /* Rotation: first around z by pole_ra, then around y' by (90-pole_dec) */
    double x_eq = -sin_ra_pole * y + cos_ra_pole * (cos_dec_pole * x + sin_dec_pole * z);
    double y_eq = cos_ra_pole * y + sin_ra_pole * (cos_dec_pole * x + sin_dec_pole * z);
    double z_eq = -sin_dec_pole * x + cos_dec_pole * z;
    
    *dec = asin(z_eq) * RAD2DEG;
    *ra = atan2(y_eq, x_eq) * RAD2DEG;
    if (*ra < 0) *ra += 360.0;
}

/**
 * @brief Generate stars for a single stream
 */
void generate_stream(FILE *fp, const StreamParams *stream, int start_id) {
    printf("Generating stream '%s' with %d stars...\n", stream->name, stream->n_stars);
    
    for (int i = 0; i < stream->n_stars; i++) {
        /* Position along stream */
        double phi1 = stream->phi1_start + 
                     rand_uniform() * (stream->phi1_end - stream->phi1_start);
        
        /* Perpendicular offset (Gaussian) */
        double phi2 = rand_gaussian() * stream->width;
        
        /* Convert to equatorial */
        double ra, dec;
        stream_to_equatorial(stream->pole_ra, stream->pole_dec, 
                            phi1, phi2, &ra, &dec);
        
        /* Proper motion along stream direction */
        /* Approximate: assume stream runs mostly E-W near equator */
        double pm_along = stream->pm_parallel + rand_gaussian() * stream->pm_dispersion;
        double pm_across = stream->pm_perp + rand_gaussian() * stream->pm_dispersion;
        
        /* Rotate to RA/Dec frame (simplified) */
        double pa = phi1 * DEG2RAD;  /* Position angle approximation */
        double pmra = pm_along * cos(pa) - pm_across * sin(pa);
        double pmdec = pm_along * sin(pa) + pm_across * cos(pa);
        
        /* Magnitude: spread around mean */
        double g_mag = stream->mean_mag + rand_gaussian() * stream->mag_spread;
        if (g_mag > 20.5) g_mag = 20.5;
        if (g_mag < 14.0) g_mag = 14.0;
        
        /* Color: narrow spread (old population) */
        double bp_rp = stream->mean_color + rand_gaussian() * stream->color_spread;
        if (bp_rp < 0.2) bp_rp = 0.2;
        if (bp_rp > 1.5) bp_rp = 1.5;
        
        fprintf(fp, "%.8f,%.8f,%.4f,%.4f,%.3f,%.3f,0.0,5000,4.0,-1.5,1,%d\n",
                ra, dec, pmra, pmdec, g_mag, bp_rp, start_id + i);
    }
}

/* ============================================================================
 * Test Scenarios
 * ============================================================================ */

/**
 * @brief Generate GD-1-like stream
 */
StreamParams create_gd1_like_stream(void) {
    StreamParams s;
    strcpy(s.name, "GD1-like");
    s.pole_ra = 34.5;
    s.pole_dec = 29.8;
    s.phi1_start = -60.0;
    s.phi1_end = 10.0;
    s.width = 0.3;
    s.pm_parallel = -12.0;
    s.pm_perp = -3.0;
    s.pm_dispersion = 0.5;
    s.mean_mag = 18.5;
    s.mag_spread = 1.0;
    s.mean_color = 0.65;
    s.color_spread = 0.1;
    s.n_stars = 500;
    return s;
}

/**
 * @brief Generate a faint, narrow stream
 */
StreamParams create_faint_stream(void) {
    StreamParams s;
    strcpy(s.name, "Faint-narrow");
    s.pole_ra = 180.0;
    s.pole_dec = 45.0;
    s.phi1_start = -30.0;
    s.phi1_end = 30.0;
    s.width = 0.2;
    s.pm_parallel = -5.0;
    s.pm_perp = 2.0;
    s.pm_dispersion = 0.3;
    s.mean_mag = 19.5;
    s.mag_spread = 0.8;
    s.mean_color = 0.7;
    s.color_spread = 0.1;
    s.n_stars = 200;
    return s;
}

/**
 * @brief Generate a bright, wide stream
 */
StreamParams create_bright_stream(void) {
    StreamParams s;
    strcpy(s.name, "Bright-wide");
    s.pole_ra = 90.0;
    s.pole_dec = -30.0;
    s.phi1_start = -20.0;
    s.phi1_end = 40.0;
    s.width = 1.0;
    s.pm_parallel = 8.0;
    s.pm_perp = -1.0;
    s.pm_dispersion = 0.8;
    s.mean_mag = 17.0;
    s.mag_spread = 1.5;
    s.mean_color = 0.75;
    s.color_spread = 0.15;
    s.n_stars = 800;
    return s;
}

/**
 * @brief Generate stream with high proper motion
 */
StreamParams create_high_pm_stream(void) {
    StreamParams s;
    strcpy(s.name, "High-PM");
    s.pole_ra = 270.0;
    s.pole_dec = 60.0;
    s.phi1_start = -15.0;
    s.phi1_end = 25.0;
    s.width = 0.4;
    s.pm_parallel = -18.0;
    s.pm_perp = 5.0;
    s.pm_dispersion = 0.4;
    s.mean_mag = 18.0;
    s.mag_spread = 1.2;
    s.mean_color = 0.6;
    s.color_spread = 0.12;
    s.n_stars = 350;
    return s;
}

/* ============================================================================
 * Main
 * ============================================================================ */

void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  -o FILE     Output file (default: synthetic_data.csv)\n");
    printf("  -n N        Number of background stars (default: 100000)\n");
    printf("  -s SEED     Random seed (default: time-based)\n");
    printf("  --scenario  Test scenario: basic, gd1, multi, stress (default: basic)\n");
    printf("  -h          Show this help\n");
}

int main(int argc, char *argv[]) {
    char output_file[256] = "synthetic_data.csv";
    int n_background = 100000;
    char scenario[32] = "basic";
    
    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            strncpy(output_file, argv[++i], sizeof(output_file) - 1);
        } else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            n_background = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            seed = atol(argv[++i]);
        } else if (strcmp(argv[i], "--scenario") == 0 && i + 1 < argc) {
            strncpy(scenario, argv[++i], sizeof(scenario) - 1);
        } else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    /* Initialize random seed if not set */
    if (seed == 12345) {
        seed = time(NULL);
    }
    
    printf("Synthetic Data Generator\n");
    printf("========================\n");
    printf("Output: %s\n", output_file);
    printf("Background stars: %d\n", n_background);
    printf("Scenario: %s\n", scenario);
    printf("Random seed: %lu\n\n", seed);
    
    FILE *fp = fopen(output_file, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s for writing\n", output_file);
        return 1;
    }
    
    /* Write header */
    fprintf(fp, "ra,dec,pmra,pmdec,g_mag,bp_rp,rv,teff,logg,feh,is_stream,target_id\n");
    
    /* Generate background */
    generate_background(fp, n_background, 20.0);
    
    int star_id = n_background;
    
    /* Generate streams based on scenario */
    if (strcmp(scenario, "basic") == 0) {
        /* Single GD-1-like stream */
        StreamParams s = create_gd1_like_stream();
        generate_stream(fp, &s, star_id);
        star_id += s.n_stars;
        
    } else if (strcmp(scenario, "gd1") == 0) {
        /* Focus on GD-1-like detection */
        StreamParams s = create_gd1_like_stream();
        s.n_stars = 1000;  /* More stars for easier detection */
        generate_stream(fp, &s, star_id);
        star_id += s.n_stars;
        
    } else if (strcmp(scenario, "multi") == 0) {
        /* Multiple streams */
        StreamParams streams[] = {
            create_gd1_like_stream(),
            create_faint_stream(),
            create_bright_stream(),
            create_high_pm_stream()
        };
        int n_streams = sizeof(streams) / sizeof(streams[0]);
        
        for (int i = 0; i < n_streams; i++) {
            generate_stream(fp, &streams[i], star_id);
            star_id += streams[i].n_stars;
        }
        
    } else if (strcmp(scenario, "stress") == 0) {
        /* Stress test: many faint streams */
        for (int i = 0; i < 10; i++) {
            StreamParams s;
            sprintf(s.name, "Stream-%d", i);
            s.pole_ra = rand_uniform() * 360.0;
            s.pole_dec = (rand_uniform() - 0.5) * 120.0;  /* Avoid poles */
            s.phi1_start = -20.0 - rand_uniform() * 20.0;
            s.phi1_end = 20.0 + rand_uniform() * 20.0;
            s.width = 0.2 + rand_uniform() * 0.5;
            s.pm_parallel = (rand_uniform() - 0.5) * 30.0;
            s.pm_perp = (rand_uniform() - 0.5) * 10.0;
            s.pm_dispersion = 0.2 + rand_uniform() * 0.5;
            s.mean_mag = 17.5 + rand_uniform() * 2.0;
            s.mag_spread = 0.8 + rand_uniform() * 0.5;
            s.mean_color = 0.5 + rand_uniform() * 0.4;
            s.color_spread = 0.1;
            s.n_stars = 100 + (int)(rand_uniform() * 300);
            
            generate_stream(fp, &s, star_id);
            star_id += s.n_stars;
        }
    }
    
    fclose(fp);
    
    printf("\nGenerated %d total stars\n", star_id);
    printf("  Background: %d\n", n_background);
    printf("  Stream: %d\n", star_id - n_background);
    printf("\nOutput written to: %s\n", output_file);
    
    /* Write ground truth file */
    char truth_file[280];
    snprintf(truth_file, sizeof(truth_file), "%s.truth", output_file);
    FILE *tfp = fopen(truth_file, "w");
    if (tfp) {
        fprintf(tfp, "# Ground truth for %s\n", output_file);
        fprintf(tfp, "# Scenario: %s\n", scenario);
        fprintf(tfp, "# Background stars: %d\n", n_background);
        fprintf(tfp, "# Stream stars: %d\n", star_id - n_background);
        fprintf(tfp, "# Stars with is_stream=1 are stream members\n");
        fclose(tfp);
        printf("Ground truth written to: %s\n", truth_file);
    }
    
    return 0;
}

