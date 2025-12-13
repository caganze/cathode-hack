# Makefile for stream detection tests (pure C)

CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c99

DATA_DIR = data
OUTPUT_DIR = output
SRC_DIR = ../src

GEN_DATA = generate_test_data
DETECTOR = $(SRC_DIR)/stream_detect

# Unit test executables
TEST_KDTREE = test_kdtree
TEST_COORDS = test_coordinates
TEST_RV = test_rv_dispersion

all: build test

# Build unit tests
unit-tests: $(TEST_KDTREE) $(TEST_COORDS) $(TEST_RV)

$(TEST_KDTREE): test_kdtree.c $(SRC_DIR)/kdtree3d.c $(SRC_DIR)/kdtree3d.h
	$(CC) $(CFLAGS) -I$(SRC_DIR) -o $@ test_kdtree.c $(SRC_DIR)/kdtree3d.c -lm

$(TEST_COORDS): test_coordinates.c $(SRC_DIR)/coordinates.o $(SRC_DIR)/utils.o
	$(CC) $(CFLAGS) -I$(SRC_DIR) -o $@ test_coordinates.c $(SRC_DIR)/coordinates.o $(SRC_DIR)/utils.o -lm

$(TEST_RV): test_rv_dispersion.c
	$(CC) $(CFLAGS) -I$(SRC_DIR) -o $@ test_rv_dispersion.c -lm

# Run unit tests
run-unit-tests: unit-tests
	@echo ""
	@echo "========================================================"
	@echo "                  UNIT TESTS"
	@echo "========================================================"
	@echo ""
	./$(TEST_KDTREE)
	@echo ""
	./$(TEST_RV)
	@echo ""
	@echo "Unit tests complete!"

$(GEN_DATA): generate_test_data.c
	$(CC) $(CFLAGS) -o $@ $< -lm

build-detector:
	$(MAKE) -C $(SRC_DIR)

build: $(GEN_DATA) build-detector

dirs:
	@mkdir -p $(DATA_DIR)
	@mkdir -p $(OUTPUT_DIR)

# Background density (stars/sq.deg) - fixed at 300
DENSITY = 300
PATCH_RADIUS = 10

# Run single test (density-based background)
define run_test
	@echo "--------------------------------------------------------"
	@echo "$(1)"
	@echo "--------------------------------------------------------"
	./$(GEN_DATA) -o $(DATA_DIR)/$(2).csv --density $(DENSITY) -r $(PATCH_RADIUS) --scenario $(3) --error $(4)
	$(DETECTOR) $(DATA_DIR)/$(2).csv -o $(OUTPUT_DIR)/$(2) -q
	@echo ""
endef

# Main test - runs ALL scenarios with ALL error levels
test: dirs build
	@echo ""
	@echo "========================================================"
	@echo "       STREAM DETECTION COMPREHENSIVE TEST SUITE"
	@echo "========================================================"
	@echo ""
	@echo "Background: $(DENSITY) stars/sq.deg, radius=$(PATCH_RADIUS) deg"
	@echo "  -> Area = $(shell echo "scale=0; 3.14159 * $(PATCH_RADIUS) * $(PATCH_RADIUS)" | bc) sq.deg"
	@echo "  -> ~$(shell echo "scale=0; $(DENSITY) * 3.14159 * $(PATCH_RADIUS) * $(PATCH_RADIUS)" | bc) background stars"
	@echo ""
	@echo "Scenarios: GD1, Diffuse (5x, 10x), Globular Cluster, Multi-stream"
	@echo "Error levels: 0%, 10%, 50%"
	@echo ""
	@# ---- GD1 Stream Tests (cold/thin) ----
	$(call run_test,GD1 Stream - No Errors,gd1_err0,gd1,0.0)
	$(call run_test,GD1 Stream - 10% Errors,gd1_err10,gd1,0.1)
	$(call run_test,GD1 Stream - 50% Errors,gd1_err50,gd1,0.5)
	@# ---- Diffuse Stream Tests (5x wider than GD1) ----
	$(call run_test,Diffuse 5x - No Errors,diffuse5_err0,diffuse5,0.0)
	$(call run_test,Diffuse 5x - 10% Errors,diffuse5_err10,diffuse5,0.1)
	$(call run_test,Diffuse 5x - 50% Errors,diffuse5_err50,diffuse5,0.5)
	@# ---- Very Diffuse Stream Tests (10x wider than GD1) ----
	$(call run_test,Diffuse 10x - No Errors,diffuse10_err0,diffuse10,0.0)
	$(call run_test,Diffuse 10x - 10% Errors,diffuse10_err10,diffuse10,0.1)
	$(call run_test,Diffuse 10x - 50% Errors,diffuse10_err50,diffuse10,0.5)
	@# ---- Globular Cluster Tests ----
	$(call run_test,Globular Cluster - No Errors,gc_err0,gc,0.0)
	$(call run_test,Globular Cluster - 10% Errors,gc_err10,gc,0.1)
	$(call run_test,Globular Cluster - 50% Errors,gc_err50,gc,0.5)
	@# ---- Multi-stream Tests ----
	$(call run_test,Multi-stream (4) - No Errors,multi_err0,multi,0.0)
	$(call run_test,Multi-stream (4) - 10% Errors,multi_err10,multi,0.1)
	$(call run_test,Multi-stream (4) - 50% Errors,multi_err50,multi,0.5)
	@# ---- Mixed (Streams + Clusters) ----
	$(call run_test,Mixed Streams+Clusters - No Errors,gc_multi_err0,gc_multi,0.0)
	$(call run_test,Mixed Streams+Clusters - 10% Errors,gc_multi_err10,gc_multi,0.1)
	$(call run_test,Mixed Streams+Clusters - 50% Errors,gc_multi_err50,gc_multi,0.5)
	@# ---- RV + Distance Tests ----
	$(call run_test,GD1 with RV+Dist - No Errors,rv_gd1_err0,rv,0.0)
	$(call run_test,GD1 with RV+Dist - 10% Errors,rv_gd1_err10,rv,0.1)
	$(call run_test,GD1 with RV+Dist - 50% Errors,rv_gd1_err50,rv,0.5)
	$(call run_test,GC with RV+Dist - No Errors,rv_gc_err0,rv_gc,0.0)
	$(call run_test,GC with RV+Dist - 10% Errors,rv_gc_err10,rv_gc,0.1)
	$(call run_test,GC with RV+Dist - 50% Errors,rv_gc_err50,rv_gc,0.5)
	$(call run_test,Multi with RV+Dist - No Errors,rv_multi_err0,rv_multi,0.0)
	$(call run_test,Multi with RV+Dist - 10% Errors,rv_multi_err10,rv_multi,0.1)
	$(call run_test,Multi with RV+Dist - 50% Errors,rv_multi_err50,rv_multi,0.5)
	@# ---- Ultra-Cold Stream Tests (RV dispersion < 5 km/s) ----
	$(call run_test,Ultra-Cold Streams - No Errors,rv_cold_err0,rv_cold,0.0)
	$(call run_test,Ultra-Cold Streams - 10% Errors,rv_cold_err10,rv_cold,0.1)
	@# ---- Mixed Cold/Hot Tests ----
	$(call run_test,Mixed Cold+Hot - No Errors,rv_mixed_err0,rv_mixed,0.0)
	$(call run_test,Mixed Cold+Hot - 10% Errors,rv_mixed_err10,rv_mixed,0.1)
	@# ---- Hot Contaminant (should be rejected) ----
	$(call run_test,Hot Contaminant - No Errors,rv_hot_err0,rv_hot,0.0)
	@# Summary
	@echo "========================================================"
	@echo "                    TEST SUMMARY"
	@echo "========================================================"
	@echo ""
	@echo "Test                           | 0% err | 10% err | 50% err"
	@echo "-------------------------------|--------|---------|--------"
	@printf "GD1 Stream (cold/thin)         |"
	@wc -l $(OUTPUT_DIR)/gd1_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/gd1_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/gd1_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Diffuse 5x                     |"
	@wc -l $(OUTPUT_DIR)/diffuse5_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/diffuse5_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/diffuse5_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Diffuse 10x                    |"
	@wc -l $(OUTPUT_DIR)/diffuse10_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/diffuse10_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/diffuse10_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Globular Cluster               |"
	@wc -l $(OUTPUT_DIR)/gc_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/gc_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/gc_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Multi-stream (4)               |"
	@wc -l $(OUTPUT_DIR)/multi_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/multi_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/multi_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Mixed (Streams+GC)             |"
	@wc -l $(OUTPUT_DIR)/gc_multi_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/gc_multi_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/gc_multi_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "GD1 + RV/Dist                  |"
	@wc -l $(OUTPUT_DIR)/rv_gd1_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_gd1_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_gd1_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "GC + RV/Dist                   |"
	@wc -l $(OUTPUT_DIR)/rv_gc_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_gc_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_gc_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Multi + RV/Dist                |"
	@wc -l $(OUTPUT_DIR)/rv_multi_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_multi_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_multi_err50/stream_candidates.csv 2>/dev/null | awk '{printf " %5d\n", $$1-1}'
	@printf "Ultra-Cold (disp<5km/s)        |"
	@wc -l $(OUTPUT_DIR)/rv_cold_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_cold_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@echo "  N/A"
	@printf "Mixed Cold+Hot                 |"
	@wc -l $(OUTPUT_DIR)/rv_mixed_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@wc -l $(OUTPUT_DIR)/rv_mixed_err10/stream_candidates.csv 2>/dev/null | awk '{printf " %5d   |", $$1-1}'
	@echo "  N/A"
	@printf "Hot Contaminant (reject!)      |"
	@wc -l $(OUTPUT_DIR)/rv_hot_err0/stream_candidates.csv 2>/dev/null | awk '{printf " %5d  |", $$1-1}'
	@echo "  N/A   |  N/A"
	@echo ""
	@echo "Results in: $(OUTPUT_DIR)/"
	@echo "========================================================"

# Quick test (smaller patch for speed)
quick: dirs build
	@echo "Quick test (density=300, radius=5 deg)..."
	./$(GEN_DATA) -o $(DATA_DIR)/quick.csv --density 300 -r 5 --scenario basic
	$(DETECTOR) $(DATA_DIR)/quick.csv -o $(OUTPUT_DIR)/quick -q

clean:
	rm -rf $(OUTPUT_DIR)/*
	rm -f $(DATA_DIR)/*.csv $(DATA_DIR)/*.txt
	rm -f $(GEN_DATA)
	rm -f $(TEST_KDTREE) $(TEST_COORDS) $(TEST_RV)

distclean: clean
	$(MAKE) -C $(SRC_DIR) clean

help:
	@echo "Stream Detection Test Suite"
	@echo ""
	@echo "  make build           - Build everything"
	@echo "  make test            - Run ALL integration tests"
	@echo "  make unit-tests      - Build unit tests"
	@echo "  make run-unit-tests  - Run unit tests (kdtree, coords, RV dispersion)"
	@echo "  make quick           - Quick test"
	@echo "  make clean           - Remove outputs"
	@echo ""
	@echo "Unit Tests:"
	@echo "  - test_kdtree       - KD-tree construction and queries"
	@echo "  - test_coordinates  - Coordinate transformations"
	@echo "  - test_rv_dispersion- RV dispersion calculations and filtering"
	@echo ""
	@echo "Integration Scenarios tested:"
	@echo "  - GD1 stream (cold/thin, width=0.3 deg, PM disp=0.5 mas/yr)"
	@echo "  - Diffuse 5x (width=1.5 deg, PM disp=2.5 mas/yr)"
	@echo "  - Diffuse 10x (width=3.0 deg, PM disp=5.0 mas/yr)"
	@echo "  - Globular Cluster (M53-like)"
	@echo "  - Multi-stream (4 streams)"
	@echo "  - Mixed (streams + clusters)"
	@echo "  - GD1 with RV/Distance"
	@echo "  - GC with RV/Distance"
	@echo "  - Multi with RV/Distance"
	@echo "  - Ultra-Cold Streams (RV disp < 5 km/s)"
	@echo "  - Mixed Cold+Hot (tests RV filtering)"
	@echo "  - Hot Contaminant (should be rejected)"
	@echo ""
	@echo "Error levels: 0%, 10%, 50%"

.PHONY: all build build-detector dirs test quick clean distclean help unit-tests run-unit-tests
