/**
 * @file coordinates.c
 * @brief Coordinate transformations and patch handling
 * 
 * Implements the coordinate transformations used in Via Machinae
 * for dividing the sky into patches and transforming to patch-local coordinates.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stream_detect.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Conversion factors */
#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

/* Galactic North Pole in equatorial coordinates (J2000) */
#define GNP_RA  192.85948    /* degrees */
#define GNP_DEC 27.12825     /* degrees */
#define GC_RA   266.405      /* Galactic center RA */
#define L_NCP   122.932      /* l of North Celestial Pole */

/**
 * @brief Convert equatorial coordinates to Galactic coordinates
 * 
 * @param ra Right ascension in degrees
 * @param dec Declination in degrees
 * @param l Output galactic longitude in degrees
 * @param b Output galactic latitude in degrees
 */
void equatorial_to_galactic(double ra, double dec, double *l, double *b) {
    double ra_rad = ra * DEG2RAD;
    double dec_rad = dec * DEG2RAD;
    double gnp_ra_rad = GNP_RA * DEG2RAD;
    double gnp_dec_rad = GNP_DEC * DEG2RAD;
    double l_ncp_rad = L_NCP * DEG2RAD;
    
    /* Calculate galactic latitude b */
    double sin_b = sin(gnp_dec_rad) * sin(dec_rad) + 
                   cos(gnp_dec_rad) * cos(dec_rad) * cos(ra_rad - gnp_ra_rad);
    *b = asin(sin_b) * RAD2DEG;
    
    /* Calculate galactic longitude l */
    double y = cos(dec_rad) * sin(ra_rad - gnp_ra_rad);
    double x = cos(gnp_dec_rad) * sin(dec_rad) - 
               sin(gnp_dec_rad) * cos(dec_rad) * cos(ra_rad - gnp_ra_rad);
    double l_rad = l_ncp_rad - atan2(y, x);
    
    /* Normalize to [0, 360) */
    *l = fmod(l_rad * RAD2DEG + 360.0, 360.0);
}

/**
 * @brief Convert Galactic coordinates to equatorial coordinates
 * 
 * @param l Galactic longitude in degrees
 * @param b Galactic latitude in degrees
 * @param ra Output right ascension in degrees
 * @param dec Output declination in degrees
 */
void galactic_to_equatorial(double l, double b, double *ra, double *dec) {
    double l_rad = l * DEG2RAD;
    double b_rad = b * DEG2RAD;
    double gnp_dec_rad = GNP_DEC * DEG2RAD;
    double l_ncp_rad = L_NCP * DEG2RAD;
    
    /* Calculate declination */
    double sin_dec = sin(gnp_dec_rad) * sin(b_rad) + 
                     cos(gnp_dec_rad) * cos(b_rad) * cos(l_ncp_rad - l_rad);
    *dec = asin(sin_dec) * RAD2DEG;
    
    /* Calculate right ascension */
    double y = cos(b_rad) * sin(l_ncp_rad - l_rad);
    double x = cos(gnp_dec_rad) * sin(b_rad) - 
               sin(gnp_dec_rad) * cos(b_rad) * cos(l_ncp_rad - l_rad);
    double ra_rad = atan2(y, x) + GNP_RA * DEG2RAD;
    
    /* Normalize to [0, 360) */
    *ra = fmod(ra_rad * RAD2DEG + 360.0, 360.0);
}

/**
 * @brief Calculate angular separation between two points on the sky
 * 
 * Uses the Haversine formula for numerical stability
 * 
 * @param ra1 RA of first point (degrees)
 * @param dec1 Dec of first point (degrees)
 * @param ra2 RA of second point (degrees)
 * @param dec2 Dec of second point (degrees)
 * @return Angular separation in degrees
 */
double angular_separation(double ra1, double dec1, double ra2, double dec2) {
    double ra1_rad = ra1 * DEG2RAD;
    double dec1_rad = dec1 * DEG2RAD;
    double ra2_rad = ra2 * DEG2RAD;
    double dec2_rad = dec2 * DEG2RAD;
    
    double dra = ra2_rad - ra1_rad;
    double ddec = dec2_rad - dec1_rad;
    
    /* Haversine formula */
    double a = sin(ddec/2) * sin(ddec/2) + 
               cos(dec1_rad) * cos(dec2_rad) * sin(dra/2) * sin(dra/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    
    return c * RAD2DEG;
}

/**
 * @brief Great circle distance (same as angular_separation)
 */
double great_circle_distance(double ra1, double dec1, double ra2, double dec2) {
    return angular_separation(ra1, dec1, ra2, dec2);
}

/**
 * @brief Transform star coordinates to patch-local coordinate system
 * 
 * Transforms (RA, Dec) to (phi, lambda) centered on the patch center.
 * This is a gnomonic (tangent plane) projection.
 * 
 * @param star Star to transform (modified in place)
 * @param center_ra Patch center RA (degrees)
 * @param center_dec Patch center Dec (degrees)
 */
void transform_to_patch_coords(Star *star, double center_ra, double center_dec) {
    double ra_rad = star->ra * DEG2RAD;
    double dec_rad = star->dec * DEG2RAD;
    double c_ra_rad = center_ra * DEG2RAD;
    double c_dec_rad = center_dec * DEG2RAD;
    
    /* Compute angular separation for cos(c) */
    double cos_c = sin(c_dec_rad) * sin(dec_rad) + 
                   cos(c_dec_rad) * cos(dec_rad) * cos(ra_rad - c_ra_rad);
    
    /* Gnomonic projection */
    double x = cos(dec_rad) * sin(ra_rad - c_ra_rad) / cos_c;
    double y = (cos(c_dec_rad) * sin(dec_rad) - 
                sin(c_dec_rad) * cos(dec_rad) * cos(ra_rad - c_ra_rad)) / cos_c;
    
    /* Store in degrees */
    star->phi = atan(x) * RAD2DEG;
    star->lambda = atan(y) * RAD2DEG;
    
    /* Compute distance from center */
    star->dist_from_center = angular_separation(star->ra, star->dec, 
                                                 center_ra, center_dec);
}

/**
 * @brief Transform proper motions to patch-local coordinate system
 * 
 * Transforms (pmra, pmdec) to (mu_phi, mu_lambda) in the patch-local frame.
 * Uses the same rotation as the coordinate transformation.
 * 
 * @param star Star to transform (modified in place)
 * @param center_ra Patch center RA (degrees)
 * @param center_dec Patch center Dec (degrees)
 */
void transform_proper_motions(Star *star, double center_ra, double center_dec) {
    double ra_rad = star->ra * DEG2RAD;
    double c_ra_rad = center_ra * DEG2RAD;
    double c_dec_rad = center_dec * DEG2RAD;
    
    /* Position angle of the great circle from star to patch center */
    double dra = c_ra_rad - ra_rad;
    double sin_pa = sin(dra) * cos(c_dec_rad);
    double cos_pa = cos(star->dec * DEG2RAD) * sin(c_dec_rad) - 
                    sin(star->dec * DEG2RAD) * cos(c_dec_rad) * cos(dra);
    double pa = atan2(sin_pa, cos_pa);
    
    /* Rotate proper motions */
    star->mu_phi = star->pmra * cos(pa) + star->pmdec * sin(pa);
    star->mu_lambda = -star->pmra * sin(pa) + star->pmdec * cos(pa);
}

/**
 * @brief Create a new patch at the specified location
 * 
 * @param center_ra Center RA (degrees)
 * @param center_dec Center Dec (degrees)
 * @param radius Patch radius (degrees)
 * @return Newly allocated Patch, or NULL on error
 */
Patch *patch_create(double center_ra, double center_dec, double radius) {
    Patch *patch = (Patch *)calloc(1, sizeof(Patch));
    if (!patch) {
        fprintf(stderr, "Error: Failed to allocate patch\n");
        return NULL;
    }
    
    patch->center_ra = center_ra;
    patch->center_dec = center_dec;
    patch->radius = radius;
    
    /* Convert to Galactic coordinates */
    equatorial_to_galactic(center_ra, center_dec, 
                          &patch->center_l, &patch->center_b);
    
    /* Allocate initial star array */
    patch->capacity = 10000;
    patch->stars = (Star *)calloc(patch->capacity, sizeof(Star));
    if (!patch->stars) {
        fprintf(stderr, "Error: Failed to allocate stars array\n");
        free(patch);
        return NULL;
    }
    
    patch->n_stars = 0;
    patch->is_valid = true;
    
    return patch;
}

/**
 * @brief Destroy a patch and free its memory
 * 
 * @param patch Patch to destroy
 */
void patch_destroy(Patch *patch) {
    if (patch) {
        if (patch->stars) {
            free(patch->stars);
        }
        free(patch);
    }
}

/**
 * @brief Add a star to a patch
 * 
 * @param patch Patch to add to
 * @param star Star to add (copied)
 * @return 0 on success, -1 on error
 */
int patch_add_star(Patch *patch, const Star *star) {
    if (!patch || !star) {
        return -1;
    }
    
    /* Resize if needed */
    if (patch->n_stars >= patch->capacity) {
        uint32_t new_capacity = patch->capacity * 2;
        Star *new_stars = (Star *)realloc(patch->stars, 
                                          new_capacity * sizeof(Star));
        if (!new_stars) {
            fprintf(stderr, "Error: Failed to reallocate stars array\n");
            return -1;
        }
        patch->stars = new_stars;
        patch->capacity = new_capacity;
    }
    
    /* Copy star and transform coordinates */
    patch->stars[patch->n_stars] = *star;
    Star *s = &patch->stars[patch->n_stars];
    s->star_idx = patch->n_stars;
    
    /* Transform to patch-local coordinates */
    transform_to_patch_coords(s, patch->center_ra, patch->center_dec);
    transform_proper_motions(s, patch->center_ra, patch->center_dec);
    
    patch->n_stars++;
    return 0;
}

/**
 * @brief Generate overlapping patches covering the sky
 * 
 * Creates patches avoiding the Galactic disk and other problematic regions.
 * Based on the Via Machinae patch layout.
 * 
 * @param patches Output array of patches (allocated by function)
 * @param n_patches Output number of patches created
 * @param cfg Configuration
 * @return 0 on success, error code otherwise
 */
int generate_patches(Patch **patches, int *n_patches, const Config *cfg) {
    if (!patches || !n_patches || !cfg) {
        return -1;
    }
    
    /* Patch spacing: patches should overlap significantly */
    double spacing = cfg->patch_radius * 0.8;  /* 80% of radius */
    
    /* Count patches first */
    int count = 0;
    double *centers_ra = NULL;
    double *centers_dec = NULL;
    int capacity = 200;
    
    centers_ra = (double *)malloc(capacity * sizeof(double));
    centers_dec = (double *)malloc(capacity * sizeof(double));
    if (!centers_ra || !centers_dec) {
        free(centers_ra);
        free(centers_dec);
        return -1;
    }
    
    /* Generate patch centers using HEALPix-like approach */
    /* For simplicity, use a regular grid in RA/Dec with spacing adjustment for Dec */
    for (double dec = -90 + spacing; dec < 90; dec += spacing) {
        /* Adjust RA spacing based on declination */
        double cos_dec = cos(dec * DEG2RAD);
        if (cos_dec < 0.01) cos_dec = 0.01;  /* Avoid division by zero near poles */
        double ra_spacing = spacing / cos_dec;
        
        for (double ra = 0; ra < 360; ra += ra_spacing) {
            /* Convert to Galactic and check if in disk */
            double l, b;
            equatorial_to_galactic(ra, dec, &l, &b);
            
            /* Skip patches too close to Galactic disk */
            if (fabs(b) < cfg->min_galactic_b) {
                continue;
            }
            
            /* Skip patches near LMC (approx RA=80, Dec=-69) */
            double lmc_dist = angular_separation(ra, dec, 80.0, -69.0);
            if (lmc_dist < 10.0) {
                continue;
            }
            
            /* Skip patches near SMC (approx RA=13, Dec=-73) */
            double smc_dist = angular_separation(ra, dec, 13.0, -73.0);
            if (smc_dist < 5.0) {
                continue;
            }
            
            /* Add patch center */
            if (count >= capacity) {
                capacity *= 2;
                centers_ra = (double *)realloc(centers_ra, capacity * sizeof(double));
                centers_dec = (double *)realloc(centers_dec, capacity * sizeof(double));
                if (!centers_ra || !centers_dec) {
                    free(centers_ra);
                    free(centers_dec);
                    return -1;
                }
            }
            
            centers_ra[count] = ra;
            centers_dec[count] = dec;
            count++;
        }
    }
    
    /* Allocate patch array */
    *patches = (Patch *)calloc(count, sizeof(Patch));
    if (!*patches) {
        free(centers_ra);
        free(centers_dec);
        return -1;
    }
    
    /* Create patches */
    for (int i = 0; i < count; i++) {
        Patch *patch = patch_create(centers_ra[i], centers_dec[i], cfg->patch_radius);
        if (!patch) {
            /* Cleanup on error */
            for (int j = 0; j < i; j++) {
                patch_destroy(&(*patches)[j]);
            }
            free(*patches);
            free(centers_ra);
            free(centers_dec);
            return -1;
        }
        (*patches)[i] = *patch;
        (*patches)[i].patch_id = i;
        free(patch);  /* Free the pointer, data was copied */
    }
    
    *n_patches = count;
    
    free(centers_ra);
    free(centers_dec);
    
    return 0;
}

/**
 * @brief Assign stars to their appropriate patches
 * 
 * Each star can belong to multiple overlapping patches.
 * 
 * @param stars Array of stars
 * @param n_stars Number of stars
 * @param patches Array of patches
 * @param n_patches Number of patches
 * @return 0 on success, error code otherwise
 */
int assign_stars_to_patches(Star *stars, uint32_t n_stars,
                           Patch **patches, int n_patches) {
    if (!stars || !patches || n_stars == 0 || n_patches == 0) {
        return -1;
    }
    
    int total_assigned = 0;
    
    for (uint32_t i = 0; i < n_stars; i++) {
        Star *star = &stars[i];
        
        for (int p = 0; p < n_patches; p++) {
            Patch *patch = patches[p];
            
            /* Check if star is within patch radius */
            double dist = angular_separation(star->ra, star->dec,
                                           patch->center_ra, patch->center_dec);
            
            if (dist <= patch->radius) {
                if (patch_add_star(patch, star) == 0) {
                    total_assigned++;
                }
            }
        }
    }
    
    printf("Assigned %d star-patch pairs\n", total_assigned);
    return 0;
}

/**
 * @brief Apply fiducial cuts to stars in a patch
 * 
 * Sets the passes_fiducial flag based on position only (no isochrone cuts).
 * Only checks that stars are within the patch radius.
 * 
 * @param patch Patch to process
 * @param cfg Configuration
 * @return Number of stars passing cuts
 */
int apply_fiducial_cuts(Patch *patch, const Config *cfg) {
    if (!patch || !cfg) {
        return -1;
    }
    
    int n_passing = 0;
    
    for (uint32_t i = 0; i < patch->n_stars; i++) {
        Star *star = &patch->stars[i];
        star->passes_fiducial = true;
        
        /* Only check angular distance from patch center */
        if (star->dist_from_center > cfg->patch_radius) {
            star->passes_fiducial = false;
            continue;
        }
        
        n_passing++;
    }
    
    return n_passing;
}

