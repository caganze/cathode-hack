/**
 * @file io.c
 * @brief File I/O routines for reading DESI data and writing results
 * 
 * Supports reading from:
 * - CSV files (simple, no dependencies)
 * - FITS files (requires cfitsio library)
 * 
 * The DESI MWS catalog can be downloaded from:
 * https://datalab.noirlab.edu/data/desi#desi-dr1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "stream_detect.h"

/* Conditional FITS support */
#ifdef USE_CFITSIO
#include <fitsio.h>
#endif

#define MAX_LINE_LEN 4096
#define MAX_FIELD_LEN 256

/* ============================================================================
 * CSV Parsing Utilities
 * ============================================================================ */

/**
 * @brief Parse a CSV line into fields
 * 
 * @param line Input line
 * @param fields Output array of field strings
 * @param max_fields Maximum number of fields
 * @return Number of fields parsed
 */
static int parse_csv_line(char *line, char **fields, int max_fields) {
    int n_fields = 0;
    char *p = line;
    
    while (*p && n_fields < max_fields) {
        /* Skip leading whitespace */
        while (*p && isspace(*p)) p++;
        
        if (*p == '\0' || *p == '\n') break;
        
        char *start = p;
        
        /* Handle quoted fields */
        if (*p == '"') {
            start = ++p;
            while (*p && *p != '"') p++;
            if (*p == '"') *p++ = '\0';
        } else {
            while (*p && *p != ',' && *p != '\n') p++;
        }
        
        if (*p == ',' || *p == '\n') {
            *p++ = '\0';
        }
        
        fields[n_fields++] = start;
    }
    
    return n_fields;
}

/**
 * @brief Find column index by name in CSV header
 * 
 * @param fields Header fields
 * @param n_fields Number of fields
 * @param name Column name to find
 * @return Column index, or -1 if not found
 */
static int find_column(char **fields, int n_fields, const char *name) {
    for (int i = 0; i < n_fields; i++) {
        /* Case-insensitive comparison */
        if (strcasecmp(fields[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

/* ============================================================================
 * CSV File Reading
 * ============================================================================ */

/**
 * @brief Read stars from a CSV file
 * 
 * Expected columns (case-insensitive):
 * Required: ra, dec, pmra, pmdec (or pm_ra, pm_dec)
 * Optional: g_mag, bp_rp, rv, teff, logg, feh, target_id
 * 
 * @param filename Path to CSV file
 * @param stars Output array of stars (allocated by function)
 * @param n_stars Output number of stars read
 * @return 0 on success, error code otherwise
 */
int read_stars_csv(const char *filename, Star **stars, uint32_t *n_stars) {
    if (!filename || !stars || !n_stars) {
        return -1;
    }
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return -1;
    }
    
    printf("Reading stars from %s...\n", filename);
    
    char line[MAX_LINE_LEN];
    char *fields[100];
    
    /* Read header */
    if (!fgets(line, MAX_LINE_LEN, fp)) {
        fprintf(stderr, "Error: Empty file\n");
        fclose(fp);
        return -1;
    }
    
    int n_header = parse_csv_line(line, fields, 100);
    
    /* Find column indices */
    int col_ra = find_column(fields, n_header, "ra");
    if (col_ra < 0) col_ra = find_column(fields, n_header, "target_ra");
    
    int col_dec = find_column(fields, n_header, "dec");
    if (col_dec < 0) col_dec = find_column(fields, n_header, "target_dec");
    
    int col_pmra = find_column(fields, n_header, "pmra");
    if (col_pmra < 0) col_pmra = find_column(fields, n_header, "pm_ra");
    if (col_pmra < 0) col_pmra = find_column(fields, n_header, "gaia_pmra");
    
    int col_pmdec = find_column(fields, n_header, "pmdec");
    if (col_pmdec < 0) col_pmdec = find_column(fields, n_header, "pm_dec");
    if (col_pmdec < 0) col_pmdec = find_column(fields, n_header, "gaia_pmdec");
    
    int col_gmag = find_column(fields, n_header, "g_mag");
    if (col_gmag < 0) col_gmag = find_column(fields, n_header, "gaia_phot_g_mean_mag");
    if (col_gmag < 0) col_gmag = find_column(fields, n_header, "phot_g_mean_mag");
    
    int col_bprp = find_column(fields, n_header, "bp_rp");
    if (col_bprp < 0) col_bprp = find_column(fields, n_header, "gaia_bp_rp");
    
    int col_rv = find_column(fields, n_header, "rv");
    if (col_rv < 0) col_rv = find_column(fields, n_header, "vrad");
    if (col_rv < 0) col_rv = find_column(fields, n_header, "radial_velocity");
    
    int col_parallax = find_column(fields, n_header, "parallax");
    if (col_parallax < 0) col_parallax = find_column(fields, n_header, "gaia_parallax");
    if (col_parallax < 0) col_parallax = find_column(fields, n_header, "plx");
    
    int col_teff = find_column(fields, n_header, "teff");
    if (col_teff < 0) col_teff = find_column(fields, n_header, "teff_gspphot");
    
    int col_teff_err = find_column(fields, n_header, "teff_err");
    
    int col_logg = find_column(fields, n_header, "logg");
    if (col_logg < 0) col_logg = find_column(fields, n_header, "logg_gspphot");
    
    int col_logg_err = find_column(fields, n_header, "logg_err");
    
    int col_feh = find_column(fields, n_header, "feh");
    if (col_feh < 0) col_feh = find_column(fields, n_header, "mh_gspphot");
    if (col_feh < 0) col_feh = find_column(fields, n_header, "fe_h");
    
    int col_feh_err = find_column(fields, n_header, "feh_err");
    if (col_feh_err < 0) col_feh_err = find_column(fields, n_header, "fe_h_err");
    
    int col_vsini = find_column(fields, n_header, "vsini");
    
    int col_distance = find_column(fields, n_header, "distance");
    if (col_distance < 0) col_distance = find_column(fields, n_header, "dist50");
    
    int col_targetid = find_column(fields, n_header, "target_id");
    if (col_targetid < 0) col_targetid = find_column(fields, n_header, "targetid");
    
    /* Check required columns */
    if (col_ra < 0 || col_dec < 0) {
        fprintf(stderr, "Error: Missing required columns (ra, dec)\n");
        fclose(fp);
        return -1;
    }
    
    printf("  Found columns: ra=%d, dec=%d, pmra=%d, pmdec=%d\n",
           col_ra, col_dec, col_pmra, col_pmdec);
    
    /* Count lines for allocation */
    uint32_t capacity = 100000;
    *stars = (Star *)calloc(capacity, sizeof(Star));
    if (!*stars) {
        fclose(fp);
        return -1;
    }
    
    *n_stars = 0;
    int lines_read = 0;
    int lines_skipped = 0;
    
    /* Read data */
    while (fgets(line, MAX_LINE_LEN, fp)) {
        lines_read++;
        
        int n_fields = parse_csv_line(line, fields, 100);
        if (n_fields <= col_ra || n_fields <= col_dec) {
            lines_skipped++;
            continue;
        }
        
        /* Resize if needed */
        if (*n_stars >= capacity) {
            capacity *= 2;
            Star *new_stars = (Star *)realloc(*stars, capacity * sizeof(Star));
            if (!new_stars) {
                fclose(fp);
                return -1;
            }
            *stars = new_stars;
        }
        
        Star *star = &(*stars)[*n_stars];
        memset(star, 0, sizeof(Star));
        
        /* Parse required fields */
        star->ra = atof(fields[col_ra]);
        star->dec = atof(fields[col_dec]);
        
        /* Skip invalid coordinates */
        if (star->ra < 0 || star->ra > 360 || star->dec < -90 || star->dec > 90) {
            lines_skipped++;
            continue;
        }
        
        /* Parse proper motions (if available) */
        if (col_pmra >= 0 && col_pmra < n_fields) {
            star->pmra = atof(fields[col_pmra]);
        }
        if (col_pmdec >= 0 && col_pmdec < n_fields) {
            star->pmdec = atof(fields[col_pmdec]);
        }
        
        /* Skip stars without proper motion */
        if (col_pmra >= 0 && col_pmdec >= 0) {
            if (star->pmra == 0 && star->pmdec == 0) {
                /* Might be valid, but often indicates missing data */
            }
        }
        
        /* Parse optional fields */
        if (col_gmag >= 0 && col_gmag < n_fields) {
            star->g_mag = atof(fields[col_gmag]);
        } else {
            star->g_mag = 20.0;  /* Default */
        }
        
        if (col_bprp >= 0 && col_bprp < n_fields) {
            star->bp_rp = atof(fields[col_bprp]);
        } else {
            star->bp_rp = 0.75;  /* Default middle of range */
        }
        
        if (col_rv >= 0 && col_rv < n_fields) {
            star->rv = atof(fields[col_rv]);
        }
        
        if (col_parallax >= 0 && col_parallax < n_fields) {
            star->parallax = atof(fields[col_parallax]);
            /* Compute distance from parallax (kpc) */
            if (star->parallax > 0.01) {
                star->distance = 1.0 / star->parallax;  /* parallax in mas -> distance in kpc */
            } else {
                star->distance = 100.0;  /* Very distant if parallax is tiny/negative */
            }
        }
        
        if (col_teff >= 0 && col_teff < n_fields) {
            star->teff = atof(fields[col_teff]);
        }
        
        if (col_teff_err >= 0 && col_teff_err < n_fields) {
            star->teff_err = atof(fields[col_teff_err]);
        }
        
        if (col_logg >= 0 && col_logg < n_fields) {
            star->logg = atof(fields[col_logg]);
        }
        
        if (col_logg_err >= 0 && col_logg_err < n_fields) {
            star->logg_err = atof(fields[col_logg_err]);
        }
        
        if (col_feh >= 0 && col_feh < n_fields) {
            star->feh = atof(fields[col_feh]);
        }
        
        if (col_feh_err >= 0 && col_feh_err < n_fields) {
            star->feh_err = atof(fields[col_feh_err]);
        }
        
        if (col_vsini >= 0 && col_vsini < n_fields) {
            star->vsini = atof(fields[col_vsini]);
        }
        
        /* Distance: prefer explicit distance column, fall back to parallax */
        if (col_distance >= 0 && col_distance < n_fields) {
            star->distance = atof(fields[col_distance]);
        } else if (col_parallax >= 0 && col_parallax < n_fields && star->parallax > 0.01) {
            star->distance = 1.0 / star->parallax;  /* parallax in mas -> distance in kpc */
        }
        
        if (col_targetid >= 0 && col_targetid < n_fields) {
            star->target_id = atoll(fields[col_targetid]);
        }
        
        star->star_idx = *n_stars;
        star->is_selected = false;
        star->passes_fiducial = true;
        
        (*n_stars)++;
        
        /* Progress update */
        if (*n_stars % 100000 == 0) {
            printf("  Read %u stars...\n", *n_stars);
        }
    }
    
    fclose(fp);
    
    printf("  Total: %u stars read, %d lines skipped\n", *n_stars, lines_skipped);
    
    return 0;
}

/* ============================================================================
 * FITS File Reading (Optional)
 * ============================================================================ */

#ifdef USE_CFITSIO

/**
 * @brief Read DESI MWS catalog from FITS file
 * 
 * @param filename Path to FITS file
 * @param stars Output array of stars (allocated by function)
 * @param n_stars Output number of stars read
 * @return 0 on success, error code otherwise
 */
int read_desi_mws_fits(const char *filename, Star **stars, uint32_t *n_stars) {
    fitsfile *fptr;
    int status = 0;
    
    printf("Reading DESI MWS catalog from %s...\n", filename);
    
    /* Open FITS file */
    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        return -1;
    }
    
    /* Move to binary table HDU */
    int hdutype;
    if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return -1;
    }
    
    /* Get number of rows */
    long nrows;
    if (fits_get_num_rows(fptr, &nrows, &status)) {
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return -1;
    }
    
    printf("  Found %ld rows\n", nrows);
    
    /* Allocate star array */
    *stars = (Star *)calloc(nrows, sizeof(Star));
    if (!*stars) {
        fits_close_file(fptr, &status);
        return -1;
    }
    
    /* Get column numbers */
    int col_ra, col_dec, col_pmra, col_pmdec;
    int col_gmag, col_bprp, col_rv, col_teff, col_logg, col_feh, col_targetid;
    
    fits_get_colnum(fptr, CASEINSEN, "TARGET_RA", &col_ra, &status);
    if (status) { status = 0; fits_get_colnum(fptr, CASEINSEN, "RA", &col_ra, &status); }
    
    fits_get_colnum(fptr, CASEINSEN, "TARGET_DEC", &col_dec, &status);
    if (status) { status = 0; fits_get_colnum(fptr, CASEINSEN, "DEC", &col_dec, &status); }
    
    fits_get_colnum(fptr, CASEINSEN, "GAIA_PMRA", &col_pmra, &status);
    if (status) { status = 0; fits_get_colnum(fptr, CASEINSEN, "PMRA", &col_pmra, &status); }
    
    fits_get_colnum(fptr, CASEINSEN, "GAIA_PMDEC", &col_pmdec, &status);
    if (status) { status = 0; fits_get_colnum(fptr, CASEINSEN, "PMDEC", &col_pmdec, &status); }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "GAIA_PHOT_G_MEAN_MAG", &col_gmag, &status);
    if (status) { status = 0; col_gmag = -1; }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "GAIA_BP_RP", &col_bprp, &status);
    if (status) { status = 0; col_bprp = -1; }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "VRAD", &col_rv, &status);
    if (status) { status = 0; col_rv = -1; }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "TEFF", &col_teff, &status);
    if (status) { status = 0; col_teff = -1; }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "LOGG", &col_logg, &status);
    if (status) { status = 0; col_logg = -1; }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "FEH", &col_feh, &status);
    if (status) { status = 0; col_feh = -1; }
    
    status = 0;
    fits_get_colnum(fptr, CASEINSEN, "TARGETID", &col_targetid, &status);
    if (status) { status = 0; col_targetid = -1; }
    
    /* Read data in chunks */
    long chunk_size = 10000;
    double *ra_buf = (double *)malloc(chunk_size * sizeof(double));
    double *dec_buf = (double *)malloc(chunk_size * sizeof(double));
    double *pmra_buf = (double *)malloc(chunk_size * sizeof(double));
    double *pmdec_buf = (double *)malloc(chunk_size * sizeof(double));
    double *gmag_buf = (double *)malloc(chunk_size * sizeof(double));
    double *bprp_buf = (double *)malloc(chunk_size * sizeof(double));
    double *rv_buf = (double *)malloc(chunk_size * sizeof(double));
    double *teff_buf = (double *)malloc(chunk_size * sizeof(double));
    double *logg_buf = (double *)malloc(chunk_size * sizeof(double));
    double *feh_buf = (double *)malloc(chunk_size * sizeof(double));
    long long *targetid_buf = (long long *)malloc(chunk_size * sizeof(long long));
    
    *n_stars = 0;
    
    for (long start = 1; start <= nrows; start += chunk_size) {
        long end = start + chunk_size - 1;
        if (end > nrows) end = nrows;
        long count = end - start + 1;
        
        /* Read columns */
        int anynul;
        double nulval = 0.0;
        
        fits_read_col(fptr, TDOUBLE, col_ra, start, 1, count, &nulval, ra_buf, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col_dec, start, 1, count, &nulval, dec_buf, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col_pmra, start, 1, count, &nulval, pmra_buf, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col_pmdec, start, 1, count, &nulval, pmdec_buf, &anynul, &status);
        
        if (col_gmag > 0) {
            fits_read_col(fptr, TDOUBLE, col_gmag, start, 1, count, &nulval, gmag_buf, &anynul, &status);
        }
        if (col_bprp > 0) {
            fits_read_col(fptr, TDOUBLE, col_bprp, start, 1, count, &nulval, bprp_buf, &anynul, &status);
        }
        if (col_rv > 0) {
            fits_read_col(fptr, TDOUBLE, col_rv, start, 1, count, &nulval, rv_buf, &anynul, &status);
        }
        if (col_teff > 0) {
            fits_read_col(fptr, TDOUBLE, col_teff, start, 1, count, &nulval, teff_buf, &anynul, &status);
        }
        if (col_logg > 0) {
            fits_read_col(fptr, TDOUBLE, col_logg, start, 1, count, &nulval, logg_buf, &anynul, &status);
        }
        if (col_feh > 0) {
            fits_read_col(fptr, TDOUBLE, col_feh, start, 1, count, &nulval, feh_buf, &anynul, &status);
        }
        if (col_targetid > 0) {
            long long nulval_ll = 0;
            fits_read_col(fptr, TLONGLONG, col_targetid, start, 1, count, &nulval_ll, targetid_buf, &anynul, &status);
        }
        
        if (status) {
            fits_report_error(stderr, status);
            break;
        }
        
        /* Copy to star array */
        for (long i = 0; i < count; i++) {
            Star *star = &(*stars)[*n_stars];
            memset(star, 0, sizeof(Star));
            
            star->ra = ra_buf[i];
            star->dec = dec_buf[i];
            star->pmra = pmra_buf[i];
            star->pmdec = pmdec_buf[i];
            
            if (col_gmag > 0) star->g_mag = gmag_buf[i];
            else star->g_mag = 20.0;
            
            if (col_bprp > 0) star->bp_rp = bprp_buf[i];
            else star->bp_rp = 0.75;
            
            if (col_rv > 0) star->rv = rv_buf[i];
            if (col_teff > 0) star->teff = teff_buf[i];
            if (col_logg > 0) star->logg = logg_buf[i];
            if (col_feh > 0) star->feh = feh_buf[i];
            if (col_targetid > 0) star->target_id = targetid_buf[i];
            
            star->star_idx = *n_stars;
            star->is_selected = false;
            star->passes_fiducial = true;
            
            (*n_stars)++;
        }
        
        if (*n_stars % 500000 == 0) {
            printf("  Read %u stars...\n", *n_stars);
        }
    }
    
    /* Cleanup */
    free(ra_buf);
    free(dec_buf);
    free(pmra_buf);
    free(pmdec_buf);
    free(gmag_buf);
    free(bprp_buf);
    free(rv_buf);
    free(teff_buf);
    free(logg_buf);
    free(feh_buf);
    free(targetid_buf);
    
    fits_close_file(fptr, &status);
    
    printf("  Total: %u stars read\n", *n_stars);
    
    return 0;
}

#else

/* Stub implementation when cfitsio not available */
int read_desi_mws_fits(const char *filename, Star **stars, uint32_t *n_stars) {
    (void)filename;
    (void)stars;
    (void)n_stars;
    fprintf(stderr, "Error: FITS support not compiled. Rebuild with -DUSE_CFITSIO and link with -lcfitsio\n");
    return -1;
}

#endif /* USE_CFITSIO */

/* ============================================================================
 * Output Writing
 * ============================================================================ */

/**
 * @brief Write stream candidates to CSV
 * 
 * @param filename Output filename
 * @param candidates Array of stream candidates
 * @param n_candidates Number of candidates
 * @return 0 on success, error code otherwise
 */
int write_stream_candidates_csv(const char *filename,
                                StreamCandidate **candidates, int n_candidates) {
    if (!candidates || n_candidates <= 0) {
        /* Create empty file */
        FILE *fp = fopen(filename, "w");
        if (fp) {
            fprintf(fp, "candidate_id,name,significance,n_protostreams,n_stars,"
                    "mean_ra,mean_dec,mean_l,mean_b,mean_pmra,mean_pmdec\n");
            fclose(fp);
        }
        printf("Wrote 0 stream candidates to %s\n", filename);
        return 0;
    }
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s for writing\n", filename);
        return -1;
    }
    
    /* Write header */
    fprintf(fp, "candidate_id,name,significance,n_protostreams,n_stars,"
            "mean_ra,mean_dec,mean_l,mean_b,mean_pmra,mean_pmdec\n");
    
    /* Count valid candidates and write */
    int n_written = 0;
    for (int i = 0; i < n_candidates; i++) {
        StreamCandidate *sc = candidates[i];
        if (!sc) continue;
        if (sc->significance <= 0) continue;
        
        fprintf(fp, "%d,%s,%.2f,%u,%u,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f\n",
                n_written + 1,
                sc->name[0] ? sc->name : "VM3-NEW",
                sc->significance,
                sc->n_protostreams,
                sc->n_stars,
                sc->mean_ra,
                sc->mean_dec,
                sc->mean_l,
                sc->mean_b,
                sc->mean_pmra,
                sc->mean_pmdec);
        n_written++;
    }
    
    fclose(fp);
    printf("Wrote %d stream candidates to %s\n", n_written, filename);
    
    return 0;
}

/**
 * @brief Write stars of a stream candidate to CSV
 * 
 * Includes all stellar properties plus anomaly detection metrics:
 * - anomaly_score: R(x) density ratio score
 * - stream_significance: overall stream significance
 * - distance_from_line: angular distance from fitted line (degrees)
 * 
 * @param filename Output filename
 * @param candidate Stream candidate
 * @return 0 on success, error code otherwise
 */
int write_stream_stars_csv(const char *filename, StreamCandidate *candidate) {
    if (!candidate || !candidate->all_stars || candidate->n_stars == 0) {
        return -1;
    }
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s for writing\n", filename);
        return -1;
    }
    
    /* Write header with full column set including anomaly metrics */
    fprintf(fp, "ra,dec,pmra,pmdec,rv,distance,g_mag,bp_rp,"
            "teff,teff_err,logg,logg_err,feh,feh_err,vsini,"
            "anomaly_score,stream_significance,dist_from_center,"
            "phi,lambda,mu_phi,mu_lambda,target_id\n");
    
    /* Write stars with all available data */
    for (uint32_t i = 0; i < candidate->n_stars; i++) {
        Star *s = candidate->all_stars[i];
        if (!s) continue;
        
        fprintf(fp, "%.8f,%.8f,%.4f,%.4f,%.2f,%.4f,%.3f,%.3f,"
                "%.1f,%.1f,%.3f,%.3f,%.3f,%.3f,%.2f,"
                "%.6f,%.2f,%.6f,"
                "%.6f,%.6f,%.4f,%.4f,%lld\n",
                s->ra, s->dec, s->pmra, s->pmdec, s->rv, s->distance,
                s->g_mag, s->bp_rp,
                s->teff, s->teff_err, s->logg, s->logg_err, 
                s->feh, s->feh_err, s->vsini,
                s->anomaly_score, candidate->significance, s->dist_from_center,
                s->phi, s->lambda, s->mu_phi, s->mu_lambda,
                (long long)s->target_id);
    }
    
    fclose(fp);
    printf("  Wrote %u stream member stars to %s\n", candidate->n_stars, filename);
    
    return 0;
}

/**
 * @brief Write patch summary
 * 
 * @param filename Output filename
 * @param patches Array of patches
 * @param n_patches Number of patches
 * @return 0 on success, error code otherwise
 */
int write_patch_summary(const char *filename, Patch **patches, int n_patches) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s for writing\n", filename);
        return -1;
    }
    
    fprintf(fp, "patch_id,center_ra,center_dec,center_l,center_b,n_stars,is_valid\n");
    
    for (int i = 0; i < n_patches; i++) {
        Patch *p = patches[i];
        fprintf(fp, "%d,%.4f,%.4f,%.4f,%.4f,%u,%d\n",
                p->patch_id, p->center_ra, p->center_dec,
                p->center_l, p->center_b, p->n_stars, p->is_valid);
    }
    
    fclose(fp);
    
    return 0;
}

