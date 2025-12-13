#!/usr/bin/env python
"""
Download the DESI DR1 Milky Way Survey (MWS) Value-Added Catalog.

The MWS VAC contains ~4 million stellar spectra with:
- Radial velocities
- Stellar parameters (Teff, logg, [Fe/H])
- Cross-matched Gaia astrometry (proper motions, parallax)

Source: https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/iron/v1.0/

Usage:
    # Download the full catalog (~13 GB FITS file)
    python download_mws.py --download -o data/
    
    # Convert to CSV for stream detection
    python download_mws.py --convert -i data/mwsall-pix-iron.fits -o data/mws_stars.csv
"""

import os
import sys
import argparse
import subprocess

# MWS VAC URLs
MWS_BASE_URL = "https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/iron/v1.0/"
MWS_CATALOG_FILE = "mwsall-pix-iron.fits"
MWS_CATALOG_URL = MWS_BASE_URL + MWS_CATALOG_FILE


def download_mws_catalog(output_dir, use_wget=True):
    """
    Download the MWS VAC FITS catalog.
    
    The file is ~13 GB so this may take a while.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, MWS_CATALOG_FILE)
    
    print("="*60)
    print("DESI DR1 MWS VAC Download")
    print("="*60)
    print(f"Source: {MWS_CATALOG_URL}")
    print(f"Output: {output_file}")
    print(f"Size: ~13 GB")
    print("")
    
    if os.path.exists(output_file):
        size_gb = os.path.getsize(output_file) / (1024**3)
        print(f"File already exists ({size_gb:.2f} GB)")
        print("Delete it first if you want to re-download.")
        return output_file
    
    if use_wget:
        print("Downloading with wget (supports resume with -c)...")
        cmd = ["wget", "-c", "--progress=bar:force", "-O", output_file, MWS_CATALOG_URL]
        try:
            subprocess.run(cmd, check=True)
            print(f"\nDownload complete: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            print(f"wget failed: {e}")
            print("Trying curl...")
    
    # Fallback to curl
    print("Downloading with curl...")
    cmd = ["curl", "-L", "-C", "-", "--progress-bar", "-o", output_file, MWS_CATALOG_URL]
    try:
        subprocess.run(cmd, check=True)
        print(f"\nDownload complete: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"curl failed: {e}")
        return None


def convert_fits_to_csv(fits_file, output_csv, columns=None, limit=None):
    """Convert MWS FITS catalog to CSV for stream detection."""
    try:
        import fitsio
        import numpy as np
    except ImportError:
        print("ERROR: fitsio and numpy required")
        print("Install with: pip install fitsio numpy")
        return False
    
    print("="*60)
    print("Converting MWS FITS to CSV")
    print("="*60)
    print(f"Input: {fits_file}")
    print(f"Output: {output_csv}")
    
    # Default columns for stream detection (positions + proper motions)
    if columns is None:
        columns = [
            ('TARGET_RA', 'ra'),
            ('TARGET_DEC', 'dec'),
            ('GAIA_PMRA', 'pmra'),
            ('GAIA_PMDEC', 'pmdec'),
            ('GAIA_PARALLAX', 'parallax'),
            ('VRAD', 'rv'),
            ('VRAD_ERR', 'rv_err'),
            ('TARGETID', 'target_id'),
        ]
    
    print(f"\nReading FITS file...")
    
    fits_cols = [c[0] for c in columns]
    try:
        data = fitsio.read(fits_file, columns=fits_cols, ext=1)
    except Exception as e:
        print(f"Error reading FITS: {e}")
        print("\nTrying to read without specifying columns...")
        data = fitsio.read(fits_file, ext=1)
    
    print(f"Total rows: {len(data):,}")
    
    available = set(data.dtype.names)
    print(f"Available columns: {len(available)}")
    
    # Filter to rows with valid proper motions
    if 'GAIA_PMRA' in available and 'GAIA_PMDEC' in available:
        valid_pm = np.isfinite(data['GAIA_PMRA']) & np.isfinite(data['GAIA_PMDEC'])
        print(f"Stars with valid proper motions: {np.sum(valid_pm):,}")
        data = data[valid_pm]
    
    if limit and len(data) > limit:
        print(f"Limiting to {limit:,} rows")
        data = data[:limit]
    
    print(f"\nWriting CSV...")
    
    with open(output_csv, 'w') as f:
        header = [csv_col for fits_col, csv_col in columns if fits_col in available]
        f.write(','.join(header) + '\n')
        
        for i, row in enumerate(data):
            vals = []
            for fits_col, csv_col in columns:
                if fits_col in available:
                    val = row[fits_col]
                    if isinstance(val, bytes):
                        val = val.decode().strip()
                    elif np.issubdtype(type(val), np.floating) and not np.isfinite(val):
                        val = ''
                    vals.append(str(val))
            
            f.write(','.join(vals) + '\n')
            
            if (i + 1) % 500000 == 0:
                print(f"  Written {i+1:,} rows...")
    
    file_size = os.path.getsize(output_csv) / (1024**2)
    print(f"\nDone! Wrote {len(data):,} stars to {output_csv} ({file_size:.1f} MB)")
    print(f"\nTo run stream detection:")
    print(f"  ./src/stream_detect {output_csv} -o output/")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Download DESI DR1 MWS Value-Added Catalog',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Source: https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/iron/v1.0/

Examples:
    python download_mws.py --download -o data/
    python download_mws.py --convert -i data/mwsall-pix-iron.fits -o data/mws_stars.csv
    python download_mws.py --download --convert -o data/
        """
    )
    
    parser.add_argument('--download', action='store_true', help='Download the MWS FITS catalog')
    parser.add_argument('--convert', action='store_true', help='Convert FITS to CSV')
    parser.add_argument('-i', '--input', type=str, help='Input FITS file (for --convert)')
    parser.add_argument('-o', '--output', type=str, default='.', help='Output directory or CSV file')
    parser.add_argument('--limit', type=int, default=None, help='Limit number of rows')
    parser.add_argument('--use-curl', action='store_true', help='Use curl instead of wget')
    
    args = parser.parse_args()
    
    if not args.download and not args.convert:
        parser.print_help()
        print("\nSpecify --download and/or --convert")
        sys.exit(1)
    
    fits_file = None
    
    if args.download:
        output_dir = args.output if not args.output.endswith('.csv') else os.path.dirname(args.output) or '.'
        fits_file = download_mws_catalog(output_dir, use_wget=not args.use_curl)
        if not fits_file:
            sys.exit(1)
    
    if args.convert:
        if args.input:
            fits_file = args.input
        elif not fits_file:
            fits_file = os.path.join(args.output, MWS_CATALOG_FILE)
        
        if not os.path.exists(fits_file):
            print(f"ERROR: FITS file not found: {fits_file}")
            print("Run with --download first, or specify -i <fits_file>")
            sys.exit(1)
        
        csv_file = args.output if args.output.endswith('.csv') else os.path.join(args.output, 'mws_stars.csv')
        
        if not convert_fits_to_csv(fits_file, csv_file, limit=args.limit):
            sys.exit(1)


if __name__ == '__main__':
    main()

