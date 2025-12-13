#!/bin/bash
###############################################################################
# Setup script for running stream detection on Sherlock
#
# This script:
# 1. Creates required directories on Oak storage
# 2. Converts the DESI LOA FITS catalog to CSV
# 3. Compiles the C code
# 4. Submits the stream detection job
#
# Usage:
#   ./setup_sherlock.sh
#
###############################################################################

set -e  # Exit on error

# Auto-detect project directory from script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Configuration - Oak storage paths
OAK_BASE="/oak/stanford/orgs/kipac/users/caganze/desi"

# FITS catalog paths
RVPIX_FILE="${OAK_BASE}/data/catalogs/rvjpix-loa.fits"
DIST_FILE="${OAK_BASE}/data/catalogs/rvsdistnn-loa-250218.fits"

# Output paths
CSV_OUTPUT="${OAK_BASE}/data/loa_stars.csv"
LOGS_DIR="${OAK_BASE}/logs"
OUTPUT_DIR="${OAK_BASE}/output"

echo "Project directory: ${PROJECT_DIR}"

echo "=============================================="
echo "Stream Detection Setup for Sherlock"
echo "=============================================="
echo ""

# Step 1: Create directories
echo "Step 1: Creating directories..."
mkdir -p "${LOGS_DIR}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OAK_BASE}/data"
echo "  Created ${LOGS_DIR}"
echo "  Created ${OUTPUT_DIR}"
echo ""

# Step 2: Check FITS files exist
echo "Step 2: Checking input FITS files..."
if [[ ! -f "${RVPIX_FILE}" ]]; then
    echo "  ERROR: RVpix file not found: ${RVPIX_FILE}"
    exit 1
fi
echo "  Found: ${RVPIX_FILE}"

if [[ ! -f "${DIST_FILE}" ]]; then
    echo "  ERROR: Distance file not found: ${DIST_FILE}"
    exit 1
fi
echo "  Found: ${DIST_FILE}"
echo ""

# Step 3: Convert FITS to CSV
echo "Step 3: Converting FITS catalog to CSV..."
if [[ -f "${CSV_OUTPUT}" ]]; then
    echo "  CSV already exists: ${CSV_OUTPUT}"
    echo "  Delete it to regenerate."
else
    echo "  Running conversion..."
    
    # Load Python module
    module load python/3.9 2>/dev/null || module load python/3.10 2>/dev/null || true
    
    # Check for required packages
    python -c "from astropy.io import fits; from astropy.table import Table" 2>/dev/null
    if [[ $? -ne 0 ]]; then
        echo "  Installing required Python packages..."
        pip install --user astropy numpy
    fi
    
    # Run conversion
    python "${PROJECT_DIR}/scripts/convert_loa_catalog.py" \
        --rvpix "${RVPIX_FILE}" \
        --dist "${DIST_FILE}" \
        --output "${CSV_OUTPUT}"
    
    if [[ $? -ne 0 ]]; then
        echo "  ERROR: Conversion failed"
        exit 1
    fi
fi

N_STARS=$(wc -l < "${CSV_OUTPUT}")
N_STARS=$((N_STARS - 1))
echo "  CSV contains ${N_STARS} stars"
echo ""

# Step 4: Compile C code
echo "Step 4: Compiling stream detection code..."
cd "${PROJECT_DIR}/src"

# Load compiler
module load gcc/12.1.0 2>/dev/null || module load gcc 2>/dev/null || true

make clean
make

if [[ ! -f "stream_detect" ]]; then
    echo "  ERROR: Compilation failed"
    exit 1
fi
echo "  Compiled: ${PROJECT_DIR}/src/stream_detect"
echo ""

# Step 5: Submit job
echo "Step 5: Submitting SLURM job..."
cd "${PROJECT_DIR}"

JOB_ID=$(sbatch --parsable scripts/run_stream_detection.sh "${CSV_OUTPUT}")
echo "  Submitted job: ${JOB_ID}"
echo ""

echo "=============================================="
echo "Setup Complete!"
echo "=============================================="
echo ""
echo "Monitor your job with:"
echo "  squeue -u ${USER}"
echo "  tail -f ${LOGS_DIR}/stream_detect_${JOB_ID}.out"
echo ""
echo "Output will be in:"
echo "  ${OUTPUT_DIR}/stream_detection_*/"

