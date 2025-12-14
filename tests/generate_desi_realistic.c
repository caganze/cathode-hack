/**
 * @file generate_desi_realistic.c
 * @brief Generate DESI-realistic synthetic data with survey artifacts
 * 
 * Creates test catalogs matching DESI DR1 characteristics:
 * - ~1.1M stars in RA [160, 280], Dec [-2, 45]
 * - Hexagonal tile/fiber footprint patterns (DESI focal plane)
 * - Realistic density variations from overlapping tiles
 * - Optional RV and distance with realistic dispersions
 * 
 * Usage:
 *   ./generate_desi_realistic --mode pos_pm         # Position + PM only
 *   ./generate_desi_realistic --mode with_dist      # + Distance
 *   ./generate_desi_realistic --mode with_rv        # + RV
 *   ./generate_desi_realistic --mode full           # All (dist + RV)
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

/* ============================================================================
 * DESI Survey Parameters
 * ============================================================================ */

/* DESI footprint for this test */
#define DESI_RA_MIN    160.0
#define DESI_RA_MAX    280.0
#define DESI_DEC_MIN   -2.0
#define DESI_DEC_MAX   45.0

/* DESI tile/fiber pattern */
#define DESI_TILE_RADIUS    1.6    /* degrees - DESI focal plane radius */
#define DESI_FIBER_RADIUS   0.05   /* degrees - individual fiber patrol radius */
#define DESI_TILE_SPACING   1.4    /* degrees - tile center spacing (overlapping) */
#define DESI_FIBERS_PER_TILE 5000  /* Approximate fibers per tile */

/* Target density */
#define TARGET_N_STARS      1100000  /* ~1.1 million stars */

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

/* ============================================================================
 * Coordinate Transformations
 * ============================================================================ */

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

/* ============================================================================
 * DESI Tile Pattern Generation
 * 
 * Creates hexagonal tile grid with realistic density variations
 * ============================================================================ */

typedef struct {
    double ra;
    double dec;
    double effective_radius;
    int n_passes;  /* Number of overlapping passes at this tile */
} DESITile;

static DESITile *g_tiles = NULL;
static int g_n_tiles = 0;

/* Generate DESI-like hexagonal tile pattern */
static void generate_tile_pattern(void) {
    /* Allocate tile array */
    int max_tiles = 10000;
    g_tiles = malloc(max_tiles * sizeof(DESITile));
    g_n_tiles = 0;
    
    /* Generate hexagonal grid of tiles */
    double dec_step = DESI_TILE_SPACING * 0.866;  /* sqrt(3)/2 for hex grid */
    int row = 0;
    
    for (double dec = DESI_DEC_MIN; dec <= DESI_DEC_MAX; dec += dec_step) {
        double cos_dec = cos(dec * DEG2RAD);
        if (cos_dec < 0.1) cos_dec = 0.1;
        double ra_step = DESI_TILE_SPACING / cos_dec;
        
        /* Offset alternating rows for hex pattern */
        double ra_offset = (row % 2 == 0) ? 0.0 : ra_step * 0.5;
        
        for (double ra = DESI_RA_MIN + ra_offset; ra <= DESI_RA_MAX; ra += ra_step) {
            if (g_n_tiles >= max_tiles) break;
            
            /* Add some randomness to tile positions (survey imperfections) */
            double tile_ra = ra + rand_gaussian() * 0.05;
            double tile_dec = dec + rand_gaussian() * 0.05;
            
            /* Clamp to footprint */
            if (tile_ra < DESI_RA_MIN) tile_ra = DESI_RA_MIN;
            if (tile_ra > DESI_RA_MAX) tile_ra = DESI_RA_MAX;
            if (tile_dec < DESI_DEC_MIN) tile_dec = DESI_DEC_MIN;
            if (tile_dec > DESI_DEC_MAX) tile_dec = DESI_DEC_MAX;
            
            /* Skip tiles in galactic plane */
            double l, b;
            equatorial_to_galactic(tile_ra, tile_dec, &l, &b);
            if (fabs(b) < 15.0) continue;
            
            g_tiles[g_n_tiles].ra = tile_ra;
            g_tiles[g_n_tiles].dec = tile_dec;
            g_tiles[g_n_tiles].effective_radius = DESI_TILE_RADIUS * (0.9 + rand_uniform() * 0.2);
            g_tiles[g_n_tiles].n_passes = 1 + (int)(rand_uniform() * 3);  /* 1-3 passes */
            g_n_tiles++;
        }
        row++;
    }
    
    printf("Generated %d DESI-like tiles\n", g_n_tiles);
}

/* Count how many tiles cover a given position (for density modulation) */
static int count_tile_coverage(double ra, double dec) {
    int count = 0;
    for (int i = 0; i < g_n_tiles; i++) {
        double dist = angular_distance(ra, dec, g_tiles[i].ra, g_tiles[i].dec);
        if (dist <= g_tiles[i].effective_radius) {
            count += g_tiles[i].n_passes;
        }
    }
    return count;
}

/* Check if position is within any tile (with edge effects) */
static int is_in_footprint(double ra, double dec, double *edge_factor) {
    *edge_factor = 1.0;
    
    for (int i = 0; i < g_n_tiles; i++) {
        double dist = angular_distance(ra, dec, g_tiles[i].ra, g_tiles[i].dec);
        if (dist <= g_tiles[i].effective_radius) {
            /* Apply edge taper (density drops near tile edge) */
            double frac = dist / g_tiles[i].effective_radius;
            if (frac > 0.8) {
                *edge_factor = 1.0 - (frac - 0.8) / 0.2 * 0.7;  /* Taper to 30% at edge */
            }
            return 1;
        }
    }
    return 0;
}

/* ============================================================================
 * Stream Definitions for Testing
 * ============================================================================ */

typedef struct {
    char name[64];
    double pole_ra, pole_dec;    /* Stream pole */
    double phi1_start, phi1_end; /* Stream extent */
    double width;                /* Angular width (degrees) */
    double pm_ra, pm_dec;        /* Mean PM (mas/yr) */
    double pm_dispersion;        /* PM dispersion (mas/yr) */
    double rv;                   /* Mean RV (km/s) */
    double rv_dispersion;        /* RV dispersion (km/s) */
    double distance;             /* Mean distance (kpc) */
    double dist_dispersion;      /* Distance dispersion (kpc) */
    int n_stars;                 /* Target number of stars */
} StreamParams;

/* Predefined test streams that pass through DESI footprint */
static StreamParams test_streams[] = {
    /* GD-1 like - cold stream, easily detectable */
    {"GD1-like", 34.5, 29.8, -40.0, 20.0, 0.3,
     -12.0, -3.0, 0.5,
     -150.0, 15.0,
     10.0, 0.8,
     3000},
    
    /* Pal5-like - ultra-cold, low RV dispersion */
    {"Pal5-like", 190.0, -30.0, -15.0, 15.0, 0.4,
     -2.1, -2.2, 0.3,
     -58.0, 2.5,  /* Very low RV dispersion! */
     23.0, 0.5,
     2000},
    
    /* Orphan-like - longer, more diffuse */
    {"Orphan-like", 310.0, 45.0, -50.0, 50.0, 1.2,
     5.0, -4.0, 1.0,
     90.0, 25.0,
     15.0, 2.0,
     4000},
    
    /* Sgr-like (leading arm piece) - very extended */
    {"Sgr-leading", 270.0, 25.0, -30.0, 60.0, 2.0,
     -2.0, -1.5, 1.5,
     -100.0, 35.0,
     25.0, 3.0,
     6000},
    
    /* Faint thin - challenging detection */
    {"Faint-thin", 220.0, 10.0, -20.0, 20.0, 0.2,
     -8.0, 2.0, 0.4,
     80.0, 20.0,
     12.0, 1.0,
     1500}
};
#define N_TEST_STREAMS (sizeof(test_streams) / sizeof(test_streams[0]))

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

/* ============================================================================
 * Data Generation
 * ============================================================================ */

typedef enum {
    MODE_POS_PM = 0,      /* Position + PM only */
    MODE_WITH_DIST = 1,   /* + Distance */
    MODE_WITH_RV = 2,     /* + RV (no distance) */
    MODE_FULL = 3         /* All (distance + RV) */
} DataMode;

static void write_header(FILE *f, DataMode mode) {
    switch (mode) {
        case MODE_POS_PM:
            fprintf(f, "ra,dec,pmra,pmdec,is_stream\n");
            break;
        case MODE_WITH_DIST:
            fprintf(f, "ra,dec,pmra,pmdec,parallax,is_stream\n");
            break;
        case MODE_WITH_RV:
            fprintf(f, "ra,dec,pmra,pmdec,rv,is_stream\n");
            break;
        case MODE_FULL:
            fprintf(f, "ra,dec,pmra,pmdec,rv,parallax,is_stream\n");
            break;
    }
}

static void write_star(FILE *f, DataMode mode, 
                       double ra, double dec, double pmra, double pmdec,
                       double rv, double distance, int is_stream) {
    double parallax = (distance > 0.1) ? 1.0 / distance : 0.0;
    
    switch (mode) {
        case MODE_POS_PM:
            fprintf(f, "%.8f,%.8f,%.4f,%.4f,%d\n", ra, dec, pmra, pmdec, is_stream);
            break;
        case MODE_WITH_DIST:
            fprintf(f, "%.8f,%.8f,%.4f,%.4f,%.6f,%d\n", 
                    ra, dec, pmra, pmdec, parallax, is_stream);
            break;
        case MODE_WITH_RV:
            fprintf(f, "%.8f,%.8f,%.4f,%.4f,%.2f,%d\n", 
                    ra, dec, pmra, pmdec, rv, is_stream);
            break;
        case MODE_FULL:
            fprintf(f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.6f,%d\n", 
                    ra, dec, pmra, pmdec, rv, parallax, is_stream);
            break;
    }
}

/* Generate background stars with DESI-like patterns */
static int generate_desi_background(FILE *f, DataMode mode, int target_n_stars) {
    printf("Generating %d DESI-like background stars...\n", target_n_stars);
    
    int generated = 0;
    int attempts = 0;
    int max_attempts = target_n_stars * 5;
    
    /* Base density (stars per attempt that get accepted) */
    double footprint_area = (DESI_RA_MAX - DESI_RA_MIN) * (DESI_DEC_MAX - DESI_DEC_MIN);
    
    while (generated < target_n_stars && attempts < max_attempts) {
        attempts++;
        
        /* Generate random position in RA/Dec box */
        double ra = DESI_RA_MIN + rand_uniform() * (DESI_RA_MAX - DESI_RA_MIN);
        double dec_sin = sin(DESI_DEC_MIN * DEG2RAD) + 
                         rand_uniform() * (sin(DESI_DEC_MAX * DEG2RAD) - sin(DESI_DEC_MIN * DEG2RAD));
        double dec = asin(dec_sin) * RAD2DEG;
        
        /* Check galactic latitude */
        double l, b;
        equatorial_to_galactic(ra, dec, &l, &b);
        if (fabs(b) < 15.0) continue;  /* Skip galactic plane */
        
        /* Check if in DESI footprint with edge effects */
        double edge_factor;
        if (!is_in_footprint(ra, dec, &edge_factor)) {
            /* Create gaps outside tiles - but keep ~20% as "serendipitous" coverage */
            if (rand_uniform() > 0.15) continue;
        }
        
        /* Modulate density by tile coverage */
        int coverage = count_tile_coverage(ra, dec);
        double accept_prob = 0.3 + coverage * 0.15;  /* More coverage = higher density */
        accept_prob *= edge_factor;
        if (rand_uniform() > accept_prob) continue;
        
        /* Generate stellar properties */
        
        /* Proper motion: smooth halo-like distribution with slight asymmetry */
        double pmra = rand_gaussian() * 4.0 + 0.3;   /* Slight mean motion */
        double pmdec = rand_gaussian() * 4.0 - 0.2;
        
        /* Add disk contamination for lower latitudes */
        if (fabs(b) < 40) {
            if (rand_uniform() < 0.3) {
                pmra += rand_gaussian() * 2.0 + 1.5 * cos(l * DEG2RAD);
                pmdec += rand_gaussian() * 1.5;
            }
        }
        
        /* RV: halo-like distribution with tails */
        double rv = rand_gaussian() * 80.0;  /* km/s */
        if (rand_uniform() < 0.1) {
            rv += (rand_uniform() > 0.5 ? 1 : -1) * 150.0;  /* High-velocity tail */
        }
        
        /* Distance: volume-weighted out to ~30 kpc */
        double dist = pow(rand_uniform(), 1.0/3.0) * 30.0;
        if (dist < 0.5) dist = 0.5;
        
        write_star(f, mode, ra, dec, pmra, pmdec, rv, dist, 0);
        generated++;
        
        if (generated % 100000 == 0) {
            printf("  Generated %d stars (%.1f%%)...\n", 
                   generated, 100.0 * generated / target_n_stars);
        }
    }
    
    printf("  Final: %d background stars\n", generated);
    return generated;
}

/* Generate stream stars */
static int generate_stream_stars(FILE *f, FILE *truth_f, DataMode mode, 
                                  StreamParams *stream) {
    printf("Generating stream '%s' with %d target stars...\n", 
           stream->name, stream->n_stars);
    
    int generated = 0;
    int attempts = 0;
    int max_attempts = stream->n_stars * 30;
    
    while (generated < stream->n_stars && attempts < max_attempts) {
        attempts++;
        
        /* Position along stream */
        double phi1 = stream->phi1_start + 
                     rand_uniform() * (stream->phi1_end - stream->phi1_start);
        double phi2 = rand_gaussian() * stream->width;
        
        double ra, dec;
        stream_to_equatorial(stream->pole_ra, stream->pole_dec, phi1, phi2, &ra, &dec);
        
        /* Check bounds */
        if (ra < DESI_RA_MIN || ra > DESI_RA_MAX) continue;
        if (dec < DESI_DEC_MIN || dec > DESI_DEC_MAX) continue;
        
        /* Check galactic latitude */
        double l, b;
        equatorial_to_galactic(ra, dec, &l, &b);
        if (fabs(b) < 15.0) continue;
        
        /* Check if in DESI footprint */
        double edge_factor;
        if (!is_in_footprint(ra, dec, &edge_factor)) {
            if (rand_uniform() > 0.15) continue;
        }
        
        /* Apply edge factor (stream stars also affected by survey incompleteness) */
        if (rand_uniform() > edge_factor) continue;
        
        /* Generate stream kinematics with tight dispersion */
        double pmra = stream->pm_ra + rand_gaussian() * stream->pm_dispersion;
        double pmdec = stream->pm_dec + rand_gaussian() * stream->pm_dispersion;
        double rv = stream->rv + rand_gaussian() * stream->rv_dispersion;
        double dist = stream->distance + rand_gaussian() * stream->dist_dispersion;
        if (dist < 0.5) dist = 0.5;
        
        write_star(f, mode, ra, dec, pmra, pmdec, rv, dist, 1);
        
        /* Write truth file entry */
        if (truth_f) {
            fprintf(truth_f, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,1,%s\n",
                    ra, dec, pmra, pmdec, rv, dist, stream->name);
        }
        
        generated++;
    }
    
    if (generated < stream->n_stars) {
        printf("  Warning: Only generated %d/%d stars (footprint constraints)\n",
               generated, stream->n_stars);
    }
    
    return generated;
}

/* ============================================================================
 * Main
 * ============================================================================ */

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Generate DESI-realistic synthetic data with survey artifacts.\n\n");
    printf("Options:\n");
    printf("  -o, --output FILE      Output CSV file (default: desi_test.csv)\n");
    printf("  -m, --mode MODE        Data mode:\n");
    printf("                           pos_pm    - Position + PM only\n");
    printf("                           with_dist - + Distance\n");
    printf("                           with_rv   - + RV\n");
    printf("                           full      - All (distance + RV)\n");
    printf("  -n, --n-stars N        Target background stars (default: 1100000)\n");
    printf("  --streams LIST         Comma-separated stream indices or 'all' (default: all)\n");
    printf("                         0=GD1, 1=Pal5, 2=Orphan, 3=Sgr, 4=Faint\n");
    printf("  --no-streams           Generate background only (null test)\n");
    printf("  -s, --seed SEED        Random seed (default: time-based)\n");
    printf("  -h, --help             Show this help\n\n");
    printf("Examples:\n");
    printf("  %s --mode pos_pm -o test_basic.csv\n", prog);
    printf("  %s --mode full -o test_full.csv --streams 0,1\n", prog);
    printf("  %s --mode with_rv --no-streams -o null_test.csv\n", prog);
}

int main(int argc, char *argv[]) {
    char output[256] = "desi_test.csv";
    DataMode mode = MODE_POS_PM;
    int n_background = TARGET_N_STARS;
    unsigned long seed = 0;
    char streams_str[256] = "all";
    int no_streams = 0;
    
    static struct option long_options[] = {
        {"output", required_argument, 0, 'o'},
        {"mode", required_argument, 0, 'm'},
        {"n-stars", required_argument, 0, 'n'},
        {"streams", required_argument, 0, 'S'},
        {"no-streams", no_argument, 0, 'N'},
        {"seed", required_argument, 0, 's'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int opt;
    while ((opt = getopt_long(argc, argv, "o:m:n:s:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'o': strncpy(output, optarg, sizeof(output)-1); break;
            case 'm':
                if (strcmp(optarg, "pos_pm") == 0) mode = MODE_POS_PM;
                else if (strcmp(optarg, "with_dist") == 0) mode = MODE_WITH_DIST;
                else if (strcmp(optarg, "with_rv") == 0) mode = MODE_WITH_RV;
                else if (strcmp(optarg, "full") == 0) mode = MODE_FULL;
                else {
                    fprintf(stderr, "Unknown mode: %s\n", optarg);
                    return 1;
                }
                break;
            case 'n': n_background = atoi(optarg); break;
            case 'S': strncpy(streams_str, optarg, sizeof(streams_str)-1); break;
            case 'N': no_streams = 1; break;
            case 's': seed = atol(optarg); break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }
    
    /* Initialize random seed */
    if (seed == 0) seed = time(NULL);
    rand_seed(seed);
    
    /* Print configuration */
    printf("DESI-Realistic Synthetic Data Generator\n");
    printf("========================================\n");
    printf("Output: %s\n", output);
    printf("Mode: %s\n", mode == MODE_POS_PM ? "pos_pm" : 
                        mode == MODE_WITH_DIST ? "with_dist" :
                        mode == MODE_WITH_RV ? "with_rv" : "full");
    printf("Target background: %d stars\n", n_background);
    printf("Footprint: RA [%.0f, %.0f], Dec [%.0f, %.0f]\n",
           DESI_RA_MIN, DESI_RA_MAX, DESI_DEC_MIN, DESI_DEC_MAX);
    printf("Streams: %s\n", no_streams ? "none" : streams_str);
    printf("Seed: %lu\n\n", seed);
    
    /* Generate tile pattern */
    generate_tile_pattern();
    
    /* Open output file */
    FILE *f = fopen(output, "w");
    if (!f) {
        fprintf(stderr, "Error: Cannot open %s for writing\n", output);
        return 1;
    }
    
    /* Open truth file */
    char truth_file[280];
    snprintf(truth_file, sizeof(truth_file), "%s", output);
    char *ext = strrchr(truth_file, '.');
    if (ext) strcpy(ext, "_truth.csv");
    else strcat(truth_file, "_truth");
    
    FILE *truth_f = fopen(truth_file, "w");
    if (truth_f) {
        fprintf(truth_f, "ra,dec,pmra,pmdec,rv,distance,is_stream,stream_name\n");
    }
    
    /* Write header */
    write_header(f, mode);
    
    /* Generate background */
    int total = generate_desi_background(f, mode, n_background);
    
    /* Generate streams */
    int stream_stars = 0;
    if (!no_streams) {
        printf("\nGenerating streams...\n");
        
        /* Parse which streams to include */
        int use_stream[N_TEST_STREAMS] = {0};
        
        if (strcmp(streams_str, "all") == 0) {
            for (size_t i = 0; i < N_TEST_STREAMS; i++) {
                use_stream[i] = 1;
            }
        } else {
            /* Parse comma-separated indices */
            char *token = strtok(streams_str, ",");
            while (token) {
                int idx = atoi(token);
                if (idx >= 0 && idx < (int)N_TEST_STREAMS) {
                    use_stream[idx] = 1;
                }
                token = strtok(NULL, ",");
            }
        }
        
        /* Generate selected streams */
        for (size_t i = 0; i < N_TEST_STREAMS; i++) {
            if (use_stream[i]) {
                int n = generate_stream_stars(f, truth_f, mode, &test_streams[i]);
                stream_stars += n;
            }
        }
    }
    
    total += stream_stars;
    
    fclose(f);
    if (truth_f) fclose(truth_f);
    free(g_tiles);
    
    printf("\n========================================\n");
    printf("Generated %d total stars\n", total);
    printf("  Background: %d\n", total - stream_stars);
    printf("  Stream: %d\n", stream_stars);
    printf("\nFiles written:\n");
    printf("  %s\n", output);
    printf("  %s\n", truth_file);
    
    return 0;
}
