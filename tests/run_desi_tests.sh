#!/bin/bash
# Run DESI-realistic stream detection tests
#
# This script tests the stream detection algorithm with synthetic data
# that mimics DESI DR1 survey characteristics including:
# - Hexagonal tile footprint patterns
# - ~1.1M stars matching your actual dataset
# - Various data modes (with/without RV, distance)
#
# Usage:
#   ./run_desi_tests.sh [quick|medium|full|null|single|all]
#
# Examples:
#   ./run_desi_tests.sh quick     # Fast test (~30s)
#   ./run_desi_tests.sh medium    # Medium tests (~5min)
#   ./run_desi_tests.sh full      # Full-scale tests (~20min)

set -e

cd "$(dirname "$0")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   DESI-Realistic Stream Detection Tests${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check/build required executables
echo -e "${YELLOW}Checking build...${NC}"
make build-detector 2>/dev/null || {
    echo -e "${RED}Error: Failed to build detector. Run 'make -C ../src' first.${NC}"
    exit 1
}

make generate_desi_realistic 2>/dev/null || {
    echo -e "${RED}Error: Failed to build data generator.${NC}"
    exit 1
}

mkdir -p data output

# Determine test mode
TEST_MODE="${1:-quick}"

case "$TEST_MODE" in
    quick)
        echo -e "${GREEN}Running quick tests (100k stars, ~30 seconds)...${NC}"
        make desi-quick
        ;;
    medium)
        echo -e "${GREEN}Running medium tests (500k stars, ~5 minutes)...${NC}"
        make desi-medium
        ;;
    full)
        echo -e "${GREEN}Running full-scale tests (1.1M stars, ~20 minutes)...${NC}"
        make desi-full
        ;;
    null)
        echo -e "${GREEN}Running null tests (false positive check)...${NC}"
        make desi-null
        ;;
    single)
        echo -e "${GREEN}Running single-stream isolation tests...${NC}"
        make desi-single
        ;;
    all)
        echo -e "${GREEN}Running ALL tests (this will take a while)...${NC}"
        make desi-all
        ;;
    compare)
        # Special mode: compare with vs without RV/distance
        echo -e "${GREEN}Running comparison: with vs without RV/distance${NC}"
        echo ""
        
        N_STARS=500000
        
        echo "Generating test data with all streams..."
        ./generate_desi_realistic -o data/compare.csv --mode full -n $N_STARS --streams all
        
        echo ""
        echo "Test 1: Position + PM only (baseline)"
        ../src/stream_detect data/compare.csv -o output/compare_pospm -q 2>&1 | tail -20
        
        echo ""
        echo "Test 2: With Distance"
        ../src/stream_detect data/compare.csv -o output/compare_dist -q --use-distance 2>&1 | tail -20
        
        echo ""
        echo "Test 3: With RV"
        ../src/stream_detect data/compare.csv -o output/compare_rv -q --use-rv 2>&1 | tail -20
        
        echo ""
        echo "Test 4: Full (Distance + RV)"
        ../src/stream_detect data/compare.csv -o output/compare_full -q --use-distance --use-rv 2>&1 | tail -20
        
        echo ""
        echo -e "${BLUE}========================================${NC}"
        echo -e "${BLUE}           COMPARISON SUMMARY${NC}"
        echo -e "${BLUE}========================================${NC}"
        echo ""
        echo "Mode          | Candidates | Notes"
        echo "--------------|------------|---------------------------"
        pospm_n=$(wc -l < output/compare_pospm/stream_candidates.csv 2>/dev/null || echo 1)
        pospm_n=$((pospm_n - 1))
        printf "pos_pm        | %10d | Baseline (position+PM)\n" $pospm_n
        
        dist_n=$(wc -l < output/compare_dist/stream_candidates.csv 2>/dev/null || echo 1)
        dist_n=$((dist_n - 1))
        printf "with_dist     | %10d | + Distance constraint\n" $dist_n
        
        rv_n=$(wc -l < output/compare_rv/stream_candidates.csv 2>/dev/null || echo 1)
        rv_n=$((rv_n - 1))
        printf "with_rv       | %10d | + RV filtering\n" $rv_n
        
        full_n=$(wc -l < output/compare_full/stream_candidates.csv 2>/dev/null || echo 1)
        full_n=$((full_n - 1))
        printf "full          | %10d | + Both (most restrictive)\n" $full_n
        
        echo ""
        echo "Expected: RV filtering should reduce false positives while"
        echo "          preserving real streams (cold RV dispersion)."
        ;;
    
    bandwidth)
        # Test different sky bandwidths to find optimal smoothing
        echo -e "${GREEN}Testing different sky bandwidth values${NC}"
        echo ""
        echo "This tests how different minimum sky bandwidths affect"
        echo "false positive rates vs stream detection sensitivity."
        echo ""
        
        N_STARS=300000
        
        echo "Generating test data..."
        ./generate_desi_realistic -o data/bw_test.csv --mode pos_pm -n $N_STARS --streams 0
        
        echo ""
        for BW in 0.5 1.0 1.5 2.0 2.5 3.0; do
            echo "Testing sky-bandwidth=${BW} deg..."
            ../src/stream_detect data/bw_test.csv -o output/bw_${BW} -q --sky-bandwidth ${BW} 2>&1 | grep -E "(candidates|Filtered)"
        done
        
        echo ""
        echo -e "${BLUE}========================================${NC}"
        echo -e "${BLUE}       SKY BANDWIDTH COMPARISON${NC}"
        echo -e "${BLUE}========================================${NC}"
        echo ""
        echo "Bandwidth | Candidates | Notes"
        echo "----------|------------|---------------------------"
        for BW in 0.5 1.0 1.5 2.0 2.5 3.0; do
            n=$(wc -l < output/bw_${BW}/stream_candidates.csv 2>/dev/null || echo 1)
            n=$((n - 1))
            note=""
            if [ "$BW" = "0.5" ]; then note="Too small - picks up tile edges"; fi
            if [ "$BW" = "2.0" ]; then note="Default - larger than tile spacing"; fi
            if [ "$BW" = "3.0" ]; then note="May over-smooth real structure"; fi
            printf "%9s | %10d | %s\n" "${BW} deg" $n "$note"
        done
        echo ""
        echo "Optimal: Bandwidth >= tile spacing (~1.4 deg) but not too large"
        ;;
    *)
        echo "Usage: $0 [quick|medium|full|null|single|all|compare|bandwidth]"
        echo ""
        echo "  quick     - Quick test (100k stars, ~30s)"
        echo "  medium    - Medium tests (500k stars, ~5min)"
        echo "  full      - Full-scale tests (1.1M stars, ~20min)"
        echo "  null      - Null test (no streams, false positive check)"
        echo "  single    - Test each stream type individually"
        echo "  all       - Run ALL tests"
        echo "  compare   - Compare with/without RV and distance"
        echo "  bandwidth - Test different sky smoothing bandwidths"
        exit 1
        ;;
esac

echo ""
echo -e "${GREEN}Tests complete!${NC}"
echo "Results are in: tests/output/"
echo ""
echo "To analyze results:"
echo "  cat output/*/stream_candidates.csv | head"
echo "  wc -l output/*/stream_candidates.csv"
