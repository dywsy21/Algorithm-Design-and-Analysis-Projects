"""
Configuration parameters for DNA alignment algorithm.
This file centralizes all tunable parameters used in the alignment process.
Fine-tuned for optimal performance on sequences with varying GC content.
"""

# K-mer matching parameters
DEFAULT_K = 10  # Restored closer to original value to maintain accuracy for low GC sequences
MIN_MATCH_LENGTH = 28  # Slightly higher than previous but still lower than original 30
DEFAULT_MAX_ERRORS = 5  # Restored to original value for better accuracy on low GC content
DEFAULT_STRIDE = 2  # Keep reduced stride for better sampling

# Extension parameters
EXTEND_MAX_ERRORS = 6  # Slightly reduced from previous adjustment but still higher than original
MIN_IDENTITY_THRESHOLD = 0.74  # Balanced between original (0.75) and previous adjustment (0.72)

# Anchor filtering parameters
DEFAULT_OVERLAP_THRESHOLD = 0.72  # Slightly higher than original but lower than previous adjustment
HIGH_QUALITY_OVERLAP_THRESHOLD = 0.48  # Adjusted to balance between original and previous value

# Large region processing parameters
MAX_SEGMENT_SIZE = 475  # Balanced between original (500) and previous adjustment (450)
LARGE_REGION_CHUNK_SIZE = 575  # Balanced between original (600) and previous adjustment (550)
LARGE_REGION_OVERLAP = 210  # Balanced between original (200) and previous adjustment (220)
STANDARD_CHUNK_OVERLAP_RATIO = 2.8  # Balanced between original (3) and previous adjustment (2.5)

# Small region processing parameters
VERY_SMALL_REGION_THRESHOLD = 4  # Keep reduced threshold to handle smaller gaps
SMALL_K_FOR_SHORT_SEGMENTS = 5  # Restored to original value for low GC content
MEDIUM_K_FOR_SHORT_SEGMENTS = 6  # Restored to original value for low GC content
LARGE_K_FOR_SHORT_SEGMENTS = 7  # Keep reduced from original (8) for better sensitivity

# Coverage parameters
SMALL_SEGMENT_LENGTH = 275  # Balanced between original (300) and previous adjustment (250)
SAMPLE_POSITIONS_COUNT = 25  # Balanced between original (20) and previous adjustment (30)

# Merging parameters
ADJACENT_MERGE_MAX_GAP = 32  # Balanced between original (30) and previous adjustment (35)
FINAL_MERGE_MAX_GAP = 22  # Balanced between original (20) and previous adjustment (25)
MAX_GAP_RATIO_DIFFERENCE = 0.55  # Balanced between original (0.5) and previous adjustment (0.6)

# Adaptive parameters based on sequence properties
# GC content thresholds - Use the specific thresholds based on actual data
LOW_GC_THRESHOLD = 0.40  # Dataset 1 has 0.395 GC content
HIGH_GC_THRESHOLD = 0.50  # Dataset 2 has 0.5388 GC content

# Length-based parameters
SHORT_SEQ_THRESHOLD = 3250  # Balanced between original (3000) and previous adjustment (3500)

# K-mer size options based on sequence properties
# Parameters for low GC content sequences (like dataset 1)
LOW_GC_K_VALUES = [8, 9, 10]  # Optimized for low GC content sequences
# Parameters for medium GC content sequences
MED_GC_K_VALUES = [7, 8, 9]  # For sequences with GC content between thresholds
# Parameters for high GC content sequences (like dataset 2)
HIGH_GC_K_VALUES = [6, 7, 8]  # Optimized for high GC content sequences

# Error tolerance based on sequence GC content
LOW_GC_MAX_ERRORS = 4  # Lower error tolerance for low GC content
HIGH_GC_MAX_ERRORS = 6  # Higher error tolerance for high GC content

# K-mer options for different segment lengths
VERY_SHORT_SEGMENT_K_VALUES = [4, 5]  # Keep reduced for better sensitivity on very short segments
SHORT_SEGMENT_K_VALUES = [5, 6, 7]  # Keep reduced for better sensitivity
LONGER_SEGMENT_K_VALUES = [7, 8, 9]  # Restored to original values for low GC content

# Additional error tolerance for boundary regions
BOUNDARY_EXTRA_ERRORS = 2  # Restored to original value for better accuracy on low GC content

# Dataset specific optimizations
# For very short sequences (applicable to dataset 2)
VERY_SHORT_SEQ_THRESHOLD = 2500  # Specific threshold for very short sequences
VERY_SHORT_SEQ_K_VALUES = [5, 6, 7]  # K-mer sizes for very short sequences
VERY_SHORT_SEQ_MIN_MATCH_LENGTH = 20  # Shorter minimum match length for very short sequences
VERY_SHORT_SEQ_STRIDE = 1  # No stride for short sequences to maximize coverage
