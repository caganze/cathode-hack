#!/usr/bin/env python
"""
Convert DESI LOA FITS catalog to CSV for the stream detection pipeline.

This script reads the DESI LOA radial velocity catalog and distance catalog,
applies quality cuts, and exports the required columns for stream detection.

Usage on Sherlock:
    module load python/3.9
    python convert_loa_catalog.py --output /oak/stanford/orgs/kipac/users/caganze/desi/data/loa_stars.csv
    
Or with custom paths:
    python convert_loa_catalog.py \
        --rvpix /oak/stanford/orgs/kipac/users/caganze/desi/data/catalogs/rvjpix-loa.fits \
        --dist /oak/stanford/orgs/kipac/users/caganze/desi/data/catalogs/rvsdistnn-loa-250218.fits \
        --output /oak/stanford/orgs/kipac/users/caganze/desi/data/loa_stars.csv
"""

import argparse
import os
import sys
import numpy as np


def read_loa_catalog(rvpix_file, dist_file, verbose=True):
    """
    Read and merge DESI LOA catalogs.
    
    Parameters
    ----------
    rvpix_file : str
        Path to rvjpix-loa.fits file
    dist_file : str
        Path to rvsdistnn-loa-*.fits distance file
    verbose : bool
        Print progress messages
        
    Returns
    -------
    table : astropy.table.Table
        Merged catalog with computed magnitudes
    """
    from astropy.io import fits
    from astropy.table import Table, hstack
    
    if verbose:
        print("="*60)
        print("Reading DESI LOA Catalog")
        print("="*60)
        print(f"RVpix file: {rvpix_file}")
        print(f"Distance file: {dist_file}")
    
    # Read rvpix catalog (may have multiple HDUs with data)
    if verbose:
        print("\nReading rvjpix-loa.fits...")
    
    with fits.open(rvpix_file) as hdul:
        if verbose:
            print(f"  Number of HDUs: {len(hdul)}")
            for i, hdu in enumerate(hdul):
                if hdu.data is not None:
                    print(f"    HDU {i}: {hdu.name}, {len(hdu.data)} rows")
        
        # Stack all data HDUs horizontally
        tables = [Table(hdu.data) for hdu in hdul if hdu.data is not None]
        stacked_table = hstack(tables)
    
    if verbose:
        print(f"  Stacked table: {len(stacked_table)} rows, {len(stacked_table.colnames)} columns")
    
    # Read distance catalog
    if verbose:
        print("\nReading distance catalog...")
    
    dt = Table.read(dist_file)
    if verbose:
        print(f"  Distance table: {len(dt)} rows, {len(dt.colnames)} columns")
    
    # Merge tables
    if verbose:
        print("\nMerging tables...")
    
    table = hstack([stacked_table, dt])
    if verbose:
        print(f"  Merged table: {len(table)} rows, {len(table.colnames)} columns")
    
    # Compute magnitudes from fluxes (based on Amanda's tutorial)
    if verbose:
        print("\nComputing magnitudes from fluxes...")
    
    # DECam magnitudes: m = 22.5 - 2.5 * log10(flux)
    for band in 'GRZ':
        flux_col = f'FLUX_{band}'
        ivar_col = f'FLUX_IVAR_{band}'
        
        if flux_col in table.colnames:
            flux = table[flux_col]
            # Handle negative/zero fluxes
            valid_flux = flux > 0
            mag = np.full(len(table), np.nan)
            mag[valid_flux] = 22.5 - 2.5 * np.log10(flux[valid_flux])
            table[f'{band.lower()}mag'] = mag
            
            # Compute magnitude errors
            if ivar_col in table.colnames:
                ivar = table[ivar_col]
                mag_err = np.full(len(table), np.nan)
                valid = valid_flux & (ivar > 0)
                mag_err[valid] = 2.5 / (np.log(10) * flux[valid] * np.sqrt(ivar[valid]))
                table[f'{band.lower()}mag_err'] = mag_err
    
    # Compute colors
    if 'gmag' in table.colnames and 'rmag' in table.colnames:
        table['g-r'] = table['gmag'] - table['rmag']
    if 'rmag' in table.colnames and 'zmag' in table.colnames:
        table['r-z'] = table['rmag'] - table['zmag']
    
    # Apply quality cuts
    if verbose:
        print("\nApplying quality cuts...")
        print(f"  Before cuts: {len(table)} stars")
    
    # Remove non-stars based on warnings and spectral type
    mask = np.ones(len(table), dtype=bool)
    
    if 'RVS_WARN' in table.colnames:
        mask &= (table['RVS_WARN'] == 0)
        if verbose:
            print(f"  After RVS_WARN==0: {np.sum(mask)} stars")
    
    if 'RR_SPECTYPE' in table.colnames:
        # Handle byte strings
        spectype = table['RR_SPECTYPE']
        if hasattr(spectype[0], 'decode'):
            spectype_str = np.array([s.decode().strip() if isinstance(s, bytes) else str(s).strip() 
                                     for s in spectype])
        else:
            spectype_str = np.array([str(s).strip() for s in spectype])
        
        mask &= (spectype_str == 'STAR')
        if verbose:
            print(f"  After RR_SPECTYPE=='STAR': {np.sum(mask)} stars")
    
    table = table[mask]
    if verbose:
        print(f"  Final catalog: {len(table)} stars")
    
    return table


def export_to_csv(table, output_file, verbose=True):
    """
    Export catalog to CSV format for the C stream detection code.
    
    Required columns: ra, dec, pmra, pmdec
    Optional columns: g_mag, bp_rp, rv, teff, logg, feh, target_id, parallax
    """
    if verbose:
        print("\n" + "="*60)
        print("Exporting to CSV")
        print("="*60)
        print(f"Output: {output_file}")
    
    # Map catalog columns to output CSV columns
    # The C code accepts various column name variants (case-insensitive)
    # Column names based on the merged LOA catalogs (rvjpix + rvsdistnn)
    column_mapping = {
        # Required: positions (from merged table, suffix _1 from first HDU)
        'ra': ['TARGET_RA_1', 'TARGET_RA', 'RA'],
        'dec': ['TARGET_DEC_1', 'TARGET_DEC', 'DEC'],
        
        # Required: proper motions (suffix _2 from second HDU or distance table)
        'pmra': ['PMRA_2', 'PMRA', 'GAIA_PMRA', 'REF_PMRA'],
        'pmdec': ['PMDEC_2', 'PMDEC', 'GAIA_PMDEC', 'REF_PMDEC'],
        
        # Optional: photometry
        'g_mag': ['GAIA_PHOT_G_MEAN_MAG', 'gmag', 'PHOT_G_MEAN_MAG'],
        'bp_rp': ['GAIA_BP_RP', 'BP_RP'],
        
        # Optional: kinematics
        'rv': ['VRAD', 'RV', 'RADIAL_VELOCITY'],
        'parallax': ['GAIA_PARALLAX', 'PARALLAX', 'REF_PARALLAX', 'PARALLAX_2'],
        
        # Optional: stellar parameters (with errors)
        'teff': ['TEFF', 'TEFF_GSPPHOT'],
        'teff_err': ['TEFF_ERR'],
        'logg': ['LOGG', 'LOGG_GSPPHOT'],
        'logg_err': ['LOGG_ERR'],
        'feh': ['FE_H', 'FEH', 'MH_GSPPHOT'],
        'feh_err': ['FE_H_ERR', 'FEH_ERR'],
        'vsini': ['VSINI'],
        
        # Optional: ID (suffix _2 for merged table)
        'target_id': ['TARGETID_2', 'TARGETID', 'TARGET_ID'],
        
        # Distance (from rvsdistnn)
        'distance': ['distance', 'dist50', 'DIST50'],
    }
    
    # Find available columns
    available_cols = set(table.colnames)
    output_cols = []
    col_sources = {}
    
    for output_name, source_options in column_mapping.items():
        for source in source_options:
            if source in available_cols:
                output_cols.append(output_name)
                col_sources[output_name] = source
                break
    
    if verbose:
        print(f"\nColumn mapping:")
        for out_name in output_cols:
            print(f"  {out_name} <- {col_sources[out_name]}")
    
    # Check required columns
    required = ['ra', 'dec', 'pmra', 'pmdec']
    missing = [c for c in required if c not in output_cols]
    if missing:
        print(f"\nERROR: Missing required columns: {missing}")
        print(f"Available columns in catalog: {sorted(available_cols)}")
        return False
    
    # Filter to valid data
    if verbose:
        print(f"\nFiltering to valid data...")
        print(f"  Initial rows: {len(table)}")
    
    valid = np.ones(len(table), dtype=bool)
    
    # Must have valid coordinates
    valid &= np.isfinite(table[col_sources['ra']])
    valid &= np.isfinite(table[col_sources['dec']])
    
    # Must have valid proper motions
    valid &= np.isfinite(table[col_sources['pmra']])
    valid &= np.isfinite(table[col_sources['pmdec']])
    
    # Reasonable coordinate ranges
    valid &= (table[col_sources['ra']] >= 0) & (table[col_sources['ra']] <= 360)
    valid &= (table[col_sources['dec']] >= -90) & (table[col_sources['dec']] <= 90)
    
    # Reasonable proper motion ranges (exclude extreme outliers)
    pm_max = 200.0  # mas/yr - reasonable for most stellar streams
    valid &= (np.abs(table[col_sources['pmra']]) < pm_max)
    valid &= (np.abs(table[col_sources['pmdec']]) < pm_max)
    
    table = table[valid]
    if verbose:
        print(f"  After filtering: {len(table)} stars")
    
    # Write CSV
    if verbose:
        print(f"\nWriting CSV file...")
    
    with open(output_file, 'w') as f:
        # Header
        f.write(','.join(output_cols) + '\n')
        
        # Data
        for i, row in enumerate(table):
            vals = []
            for out_name in output_cols:
                source = col_sources[out_name]
                val = row[source]
                
                # Handle different data types
                if isinstance(val, bytes):
                    val = val.decode().strip()
                elif np.issubdtype(type(val), np.floating):
                    if not np.isfinite(val):
                        val = ''
                    else:
                        # Use appropriate precision
                        if out_name in ['ra', 'dec']:
                            val = f'{val:.8f}'
                        elif out_name in ['pmra', 'pmdec', 'parallax']:
                            val = f'{val:.4f}'
                        elif out_name in ['g_mag', 'bp_rp', 'rv', 'logg', 'feh']:
                            val = f'{val:.3f}'
                        elif out_name in ['teff']:
                            val = f'{val:.1f}'
                        elif out_name in ['distance']:
                            val = f'{val:.4f}'
                        else:
                            val = str(val)
                else:
                    val = str(val)
                
                vals.append(val)
            
            f.write(','.join(vals) + '\n')
            
            if verbose and (i + 1) % 500000 == 0:
                print(f"  Written {i+1:,} rows...")
    
    file_size_mb = os.path.getsize(output_file) / (1024**2)
    if verbose:
        print(f"\nDone! Wrote {len(table):,} stars to {output_file}")
        print(f"File size: {file_size_mb:.1f} MB")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Convert DESI LOA FITS catalog to CSV for stream detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Using default Oak paths
    python convert_loa_catalog.py --output loa_stars.csv
    
    # Custom paths
    python convert_loa_catalog.py \\
        --rvpix /path/to/rvjpix-loa.fits \\
        --dist /path/to/rvsdistnn-loa.fits \\
        --output loa_stars.csv

After conversion, run stream detection with:
    ./stream_detect loa_stars.csv -o output/
        """
    )
    
    # Default paths on Oak storage
    default_oak = '/oak/stanford/orgs/kipac/users/caganze/desi/data/catalogs'
    
    parser.add_argument('--rvpix', type=str,
                       default=f'{default_oak}/rvjpix-loa.fits',
                       help='Path to rvjpix-loa.fits file')
    parser.add_argument('--dist', type=str,
                       default=f'{default_oak}/rvsdistnn-loa-250218.fits',
                       help='Path to rvsdistnn distance catalog')
    parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output CSV file path')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Reduce output verbosity')
    
    args = parser.parse_args()
    
    verbose = not args.quiet
    
    # Check input files exist
    for path, name in [(args.rvpix, 'rvpix'), (args.dist, 'distance')]:
        if not os.path.exists(path):
            print(f"ERROR: {name} file not found: {path}")
            sys.exit(1)
    
    # Read and merge catalogs
    try:
        table = read_loa_catalog(args.rvpix, args.dist, verbose=verbose)
    except Exception as e:
        print(f"ERROR reading catalog: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Export to CSV
    if not export_to_csv(table, args.output, verbose=verbose):
        sys.exit(1)
    
    print("\n" + "="*60)
    print("Conversion complete!")
    print("="*60)
    print(f"\nTo run stream detection on Sherlock:")
    print(f"  sbatch scripts/run_stream_detection.sh {args.output}")


if __name__ == '__main__':
    main()

