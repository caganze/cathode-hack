# Via Machinae Stream Detector

A C implementation of stellar stream detection using anomaly detection methods, inspired by the Via Machinae algorithm described in [arXiv:2509.08064v1](https://arxiv.org/abs/2509.08064) "Via Machinae 3.0: A search for stellar streams in Gaia with the CATHODE algorithm".

## Overview

This tool searches for stellar streams in the DESI DR1 Milky Way Survey (MWS) catalog using the Via Machinae methodology:

1. **Sky Division**: Divides the sky into overlapping circular patches
2. **Anomaly Detection**: Uses density estimation to find overdense regions (anomalies) in proper motion and position space
3. **Line Detection**: Applies Hough transform to find line-like structures among anomalous stars
4. **Clustering**: Groups detections into protoclusters, protostreams, and final stream candidates

## Building

### Basic Build

```bash
cd src
make
```

### With FITS Support

To read FITS files directly, you need the cfitsio library:

```bash
# Install cfitsio (macOS)
brew install cfitsio

# Or on Ubuntu/Debian
sudo apt-get install libcfitsio-dev

# Edit Makefile to uncomment FITS support lines, then:
make clean
make
```

## Getting DESI DR1 MWS Data

### Prerequisites

1. **DESI Credentials**: Create `~/.desi_http_user` with your credentials:
   ```
   username:password
   ```
   Get credentials at: https://data.desi.lbl.gov/doc/access/

2. **Install desida** (official DESI data tool from https://github.com/desihub/desida):
   ```bash
   pip install git+https://github.com/desihub/desida.git
   ```

### Option 1: Quick Setup Script (Recommended)

```bash
# Download a test subset (~1000 stars around healpix 23040)
chmod +x scripts/setup_desi_data.sh
./scripts/setup_desi_data.sh data/desi/

# Run stream detection
./src/stream_detect data/desi/mws_stars.csv -o output/
```

### Option 2: Using desi_get_dr_subset directly

```bash
# Download DR1 data subset
desi_get_dr_subset --dr dr1 --ra 56 --dec -9 --radius 0.5 --no-tiles

# Extract MWS stars to CSV
python scripts/download_desi_mws.py \
    --extract-only \
    --input tiny_dr1/ \
    --output mws_stars.csv
```

### Option 3: Full MWS Catalog via Data Lab

For the full ~4 million star catalog:

```bash
# Install Data Lab client
pip install astro-datalab

# Download full catalog (may take several minutes)
python scripts/download_desi_mws.py --full-catalog -o data/desi_mws_full.csv
```

Or query directly at https://datalab.noirlab.edu/query.php:

```sql
SELECT target_ra as ra, target_dec as dec, 
       gaia_pmra as pmra, gaia_pmdec as pmdec,
       gaia_phot_g_mean_mag as g_mag,
       gaia_phot_bp_mean_mag - gaia_phot_rp_mean_mag as bp_rp
FROM desi_dr1.zpix_redshifts
WHERE spectype = 'STAR' AND gaia_pmra IS NOT NULL
```

### Option 4: Synthetic Test Data

Generate synthetic data with injected streams for testing:

```bash
cd tests
python3 generate_test_data_simple.py -o data/test.csv --scenario multi -n 50000
```

## Usage

### Basic Usage

```bash
./stream_detect input_stars.csv -o output_dir
```

### Command Line Options

```
Usage: stream_detect [options] <input_file>

Options:
  -h, --help           Show help message
  -v, --version        Show version
  -o, --output DIR     Output directory (default: ./output)
  -p, --patch-radius R Patch radius in degrees (default: 15.0)
  -b, --min-b B        Minimum galactic |b| (default: 20.0)
  -s, --sig-cut S      Protocluster significance cut (default: 8.8)
  -k, --knn K          Use KNN density with K neighbors (default: KDE)
  --histogram          Use histogram density (faster, less accurate)
  -q, --quiet          Quiet mode
```

### Input Format

CSV format with columns:
- Required: `ra`, `dec` (or `target_ra`, `target_dec`)
- Recommended: `pmra`, `pmdec` (or `gaia_pmra`, `gaia_pmdec`)
- Optional: `g_mag`, `bp_rp`, `rv`, `teff`, `logg`, `feh`, `target_id`

### Output

The tool generates:
- `stream_candidates.csv`: Summary of all detected stream candidates
- `stream_XXX_stars.csv`: Stars belonging to each stream candidate

## Algorithm Details

### 1. Patch Division

The sky is divided into overlapping circular patches of 15° radius, excluding:
- The Galactic disk (|b| < 20°)
- The Large and Small Magellanic Clouds

### 2. Signal Regions

Within each patch, proper motion windows of 6 mas/yr width define signal regions (SRs). The complementary regions serve as sidebands (SBs) for background estimation.

### 3. Anomaly Scores

Stars in the SR receive anomaly scores R(x) = p_data(x) / p_background(x) using density estimation:
- **KDE** (default): Kernel Density Estimation with Silverman bandwidth
- **KNN**: K-Nearest Neighbors density estimation
- **Histogram**: Fast but coarse density estimation

### 4. Line Detection

The top 100 most anomalous stars in each Region of Interest (ROI) are analyzed using the Hough transform to find line-like structures.

### 5. Clustering Hierarchy

1. **ROIs → Protoclusters**: ROIs with compatible lines and proper motions are grouped
2. **Protoclusters → Protostreams**: Duplicate detection removes redundant protoclusters
3. **Protostreams → Stream Candidates**: Merging across patches

### 6. Quality Cuts

- Significance cut (σ > 8.8)
- Fraction of dim stars cut (f_dim > 0.25 with g_c = 19.1)
- Edge protostream removal

## References

- Hallin et al. (2024), "Via Machinae 3.0: A search for stellar streams in Gaia with the CATHODE algorithm", [arXiv:2509.08064](https://arxiv.org/abs/2509.08064)
- DESI Collaboration (2024), "DESI Data Release 1"
- Original Via Machinae: Shih et al. (2021, 2023)

## License

This code is provided for research and educational purposes. Please cite the original Via Machinae papers if you use this code in your research.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

