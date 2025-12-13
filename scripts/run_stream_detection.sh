#!/bin/bash
#SBATCH --job-name=stream_detect
#SBATCH --partition=kipac
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/oak/stanford/orgs/kipac/users/caganze/desi/logs/stream_detect_%j.out
#SBATCH --error=/oak/stanford/orgs/kipac/users/caganze/desi/logs/stream_detect_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=caganze@stanford.edu

###############################################################################
# DESI LOA Stream Detection Pipeline on Sherlock
#
# This script runs the Via Machinae-inspired stream detection algorithm
# on DESI LOA catalog data.
#
# Usage:
#   sbatch run_stream_detection.sh /path/to/input.csv
#   sbatch run_stream_detection.sh  # uses default path
#
# With sky region and distance cuts:
#   sbatch --export=RA_MIN=100,RA_MAX=200,DEC_MIN=-30,DEC_MAX=30 run_stream_detection.sh
#   sbatch --export=DIST_MIN=5,DIST_MAX=50 run_stream_detection.sh
#   sbatch --export=RA_MIN=100,RA_MAX=200,DIST_MIN=5,DIST_MAX=50,SIG_CUT=10 run_stream_detection.sh
#
###############################################################################

# Load required modules
module load gcc/12.1.0

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OAK_DIR="/oak/stanford/orgs/kipac/users/caganze/desi"

# Input file (can be passed as argument or use default)
INPUT_FILE="${1:-${OAK_DIR}/data/loa_stars.csv}"

# Output directory
OUTPUT_DIR="${OAK_DIR}/output/stream_detection_$(date +%Y%m%d_%H%M%S)"

# Create directories
mkdir -p "${OAK_DIR}/logs"
mkdir -p "${OUTPUT_DIR}"

echo "=============================================="
echo "DESI LOA Stream Detection Pipeline"
echo "=============================================="
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo ""
echo "Input: ${INPUT_FILE}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Check input file exists
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    echo ""
    echo "Did you run the conversion script first?"
    echo "  python scripts/convert_loa_catalog.py --output ${INPUT_FILE}"
    exit 1
fi

# Build the stream detection code if needed
cd "${PROJECT_DIR}/src"

if [[ ! -f "stream_detect" ]] || [[ "stream_detect" -ot "main.c" ]]; then
    echo "Building stream detection code..."
    make clean
    make
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Build failed"
        exit 1
    fi
    echo "Build complete."
    echo ""
fi

# Count input stars
N_STARS=$(wc -l < "${INPUT_FILE}")
N_STARS=$((N_STARS - 1))  # Subtract header
echo "Input catalog contains ${N_STARS} stars"
echo ""

# Run the pipeline
echo "Starting stream detection..."
echo "=============================================="

# Optional parameters (can be set via sbatch --export or uncomment below)
# RA_MIN=100      # Minimum RA (degrees)
# RA_MAX=200      # Maximum RA (degrees)
# DEC_MIN=-30     # Minimum Dec (degrees)
# DEC_MAX=30      # Maximum Dec (degrees)
# DIST_MIN=5      # Minimum distance (kpc)
# DIST_MAX=50     # Maximum distance (kpc)
# SIG_CUT=8.0     # Significance cutoff

# Use provided SIG_CUT or default to 8.0
SIG_CUT="${SIG_CUT:-8.0}"

# Build the command
CMD="./stream_detect ${INPUT_FILE} --output ${OUTPUT_DIR} --patch-radius 1.78 --min-b 20.0 --sig-cut ${SIG_CUT} --knn 15"

# Add optional cuts if defined
[[ -n "${RA_MIN:-}" ]] && CMD+=" --ra-min ${RA_MIN}"
[[ -n "${RA_MAX:-}" ]] && CMD+=" --ra-max ${RA_MAX}"
[[ -n "${DEC_MIN:-}" ]] && CMD+=" --dec-min ${DEC_MIN}"
[[ -n "${DEC_MAX:-}" ]] && CMD+=" --dec-max ${DEC_MAX}"
[[ -n "${DIST_MIN:-}" ]] && CMD+=" --dist-min ${DIST_MIN}"
[[ -n "${DIST_MAX:-}" ]] && CMD+=" --dist-max ${DIST_MAX}"

echo "Running: ${CMD}"
time ${CMD}

# Note: sig-cut=8.0 is aggressive for dense DESI backgrounds
# Adjust if too few/many candidates are found:
#   - Lower (5.0-7.0) for more candidates (higher false positive rate)
#   - Higher (9.0-12.0) for fewer, more reliable candidates
#
# Sky region options:
#   --ra-min, --ra-max    : RA range in degrees (0-360)
#   --dec-min, --dec-max  : Dec range in degrees (-90 to 90)
#   --dist-min, --dist-max: Distance range in kpc

EXIT_CODE=$?

echo ""
echo "=============================================="
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "Pipeline completed successfully!"
    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}/"
    echo ""
    echo "Stream candidates:"
    if [[ -f "${OUTPUT_DIR}/stream_candidates.csv" ]]; then
        head -20 "${OUTPUT_DIR}/stream_candidates.csv"
    fi
else
    echo "Pipeline failed with exit code ${EXIT_CODE}"
fi

echo ""
echo "Finished: $(date)"

