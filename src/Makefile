# Makefile for Via Machinae Stream Detector
# 
# Stellar stream detection for DESI DR1 MWS data
# Based on arXiv:2509.08064v1

CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c99
LDFLAGS = -lm

# Optional FITS support (uncomment to enable)
# CFLAGS += -DUSE_CFITSIO
# LDFLAGS += -lcfitsio

# Debug build (uncomment for debugging)
# CFLAGS = -Wall -Wextra -g -O0 -std=c99 -DDEBUG

# Source files
SRCS = main.c coordinates.c density_estimator.c hough_transform.c clustering.c io.c utils.c kdtree3d.c
OBJS = $(SRCS:.c=.o)

# Output binary
TARGET = stream_detect

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Compile
%.o: %.c stream_detect.h
	$(CC) $(CFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(TARGET)

# Install (to /usr/local/bin)
install: $(TARGET)
	install -m 755 $(TARGET) /usr/local/bin/

# Uninstall
uninstall:
	rm -f /usr/local/bin/$(TARGET)

# Run tests with sample data
test: $(TARGET)
	@echo "Running basic test..."
	@if [ -f ../data/test_stars.csv ]; then \
		./$(TARGET) ../data/test_stars.csv -o ../output/test; \
	else \
		echo "No test data found. Please download DESI MWS data first."; \
	fi

# Create test data (small sample)
sample_data:
	@echo "Creating sample test data..."
	@mkdir -p ../data
	@echo "ra,dec,pmra,pmdec,g_mag,bp_rp" > ../data/test_stars.csv
	@awk 'BEGIN { \
		srand(); \
		for (i = 0; i < 10000; i++) { \
			ra = rand() * 360; \
			dec = (rand() - 0.5) * 180; \
			pmra = (rand() - 0.5) * 20; \
			pmdec = (rand() - 0.5) * 20; \
			g = 15 + rand() * 5; \
			bp_rp = 0.3 + rand() * 1.0; \
			printf "%.6f,%.6f,%.4f,%.4f,%.3f,%.3f\n", ra, dec, pmra, pmdec, g, bp_rp; \
		} \
	}' >> ../data/test_stars.csv
	@echo "Created ../data/test_stars.csv with 10000 random stars"

# Help
help:
	@echo "Via Machinae Stream Detector - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build the stream detector (default)"
	@echo "  clean        - Remove build artifacts"
	@echo "  install      - Install to /usr/local/bin"
	@echo "  uninstall    - Remove from /usr/local/bin"
	@echo "  test         - Run tests with sample data"
	@echo "  sample_data  - Generate random test data"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "To enable FITS support, edit Makefile and uncomment:"
	@echo "  CFLAGS += -DUSE_CFITSIO"
	@echo "  LDFLAGS += -lcfitsio"
	@echo ""
	@echo "For debug build, uncomment the debug CFLAGS line"

.PHONY: all clean install uninstall test sample_data help

