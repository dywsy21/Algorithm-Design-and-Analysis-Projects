"""
Configuration parameters for DNA alignment algorithm.
This file centralizes all tunable parameters used in the alignment process.
"""

# K-mer matching parameters
DEFAULT_K = 11  # Default k-mer size for exact matching
MIN_MATCH_LENGTH = 30  # Minimum length for matches to be considered
DEFAULT_MAX_ERRORS = 5  # Default maximum number of errors allowed in matches
DEFAULT_STRIDE = 3  # Default stride for sampling positions

# Extension parameters
EXTEND_MAX_ERRORS = 5  # Maximum errors allowed when extending matches
MIN_IDENTITY_THRESHOLD = 0.75  # Minimum identity threshold for extended matches

# Anchor filtering parameters
DEFAULT_OVERLAP_THRESHOLD = 0.7  # Default overlap threshold for filtering anchors
HIGH_QUALITY_OVERLAP_THRESHOLD = 0.5  # Lower overlap threshold for higher quality anchors

# Large region processing parameters
MAX_SEGMENT_SIZE = 500  # Maximum size of segments for divide-and-conquer
LARGE_REGION_CHUNK_SIZE = 600  # Chunk size for large regions
LARGE_REGION_OVERLAP = 200  # Overlap size for large region chunks
STANDARD_CHUNK_OVERLAP_RATIO = 3  # 1/3 overlap for standard chunks

# Small region processing parameters
VERY_SMALL_REGION_THRESHOLD = 5  # Threshold for very small regions
SMALL_K_FOR_SHORT_SEGMENTS = 5  # K-mer size for very short segments
MEDIUM_K_FOR_SHORT_SEGMENTS = 6  # K-mer size for medium-short segments
LARGE_K_FOR_SHORT_SEGMENTS = 8  # K-mer size for larger segments

# Coverage parameters
SMALL_SEGMENT_LENGTH = 300  # Size of segments for uncovered regions
SAMPLE_POSITIONS_COUNT = 20  # Number of positions to sample in reference

# Merging parameters
ADJACENT_MERGE_MAX_GAP = 30  # Maximum gap for merging adjacent segments
FINAL_MERGE_MAX_GAP = 20  # Maximum gap for final merging
MAX_GAP_RATIO_DIFFERENCE = 0.5  # Maximum allowed ratio difference between query and ref gaps

# Adaptive parameters based on sequence properties
# GC content thresholds
HIGH_GC_CONTENT = 0.55  # Threshold for high GC content

# Length-based parameters
SHORT_SEQ_THRESHOLD = 3000  # Threshold for short sequences

# K-mer size options based on sequence properties
SHORT_SEQ_LOW_GC_K_VALUES = [6, 7, 8]  # K values for short sequences with low GC
SHORT_SEQ_HIGH_GC_K_VALUES = [7, 8, 9]  # K values for short sequences with high GC
LONG_SEQ_LOW_GC_K_VALUES = [8, 9, 10]  # K values for long sequences with low GC
LONG_SEQ_HIGH_GC_K_VALUES = [9, 10, 11]  # K values for long sequences with high GC

# Error tolerance based on sequence length
SHORT_SEQ_MAX_ERRORS = 3  # Max errors for short sequences
LONG_SEQ_MAX_ERRORS = 5  # Max errors for long sequences

# K-mer options for different segment lengths
VERY_SHORT_SEGMENT_K_VALUES = [5, 6]  # K values for very short segments (<100)
SHORT_SEGMENT_K_VALUES = [6, 7, 8]  # K values for short segments (<300)
LONGER_SEGMENT_K_VALUES = [7, 8, 9]  # K values for longer segments (>=300)

# Additional error tolerance for boundary regions
BOUNDARY_EXTRA_ERRORS = 2  # Extra errors allowed for boundary regions
