# Core parameters
KMER_SIZE = 15
MIN_SEGMENT_LENGTH = 30  # Required by scoring function
MAX_EDIT_RATE = 0.1      # Required by scoring function

# Multi-strategy parameters
SMALL_KMER_SIZE = 8      # For sensitive detection
LARGE_KMER_SIZE = 18     # Large k-mer for high specificity
DENSE_SHIFT = 1          # Dense sampling for better coverage
SPARSE_SHIFT = 3         # Sparse sampling for speed

# Sequence length threshold for strategy selection
SHORT_SEQUENCE_THRESHOLD = 3000

# Short sequence k-mer sampling parameters
T2_SMALL_KMER_SHIFT = 3  # Shift for very small k-mers (k<=5)
T2_MEDIUM_KMER_SHIFT = 1 # Shift for medium k-mers (k<=8)

# Long sequence boundary extension parameters
T1_BOUNDARY_TOLERANCE = 10    # Boundary proximity tolerance for short sequences
T1_RC_MISMATCH_RATE = 0.15    # Allowed mismatch rate for RC extension in short sequences
T1_RC_MIN_LENGTH = 10         # Minimum length before applying mismatch rate

# Long sequence gap filling parameters
T1_LARGE_GAP_MAX = 60         # Maximum large gap size for long sequences
T1_LARGE_GAP_MATCH_RATE = 0.85 # Required match rate for large gaps
T1_MERGED_SEGMENT_EDIT_RATE = 0.12 # More lenient edit rate for merged segments

# Small gap filling parameters for short sequences
T2_SMALL_GAP_MAX = 20         # Maximum small gap size for short sequences
T2_SMALL_GAP_MISMATCH_RATE = 0.15 # Allowed mismatch rate for small gaps

# T2 parameters (short sequences) - current best performing values
T2_KMER_SIZES = [6, 7, 8, 9, 11, 13, 15]  # K-mer sizes for short sequences
T2_EXTENSION_LIMIT_SHORT = 74  # Extension limit for short sequences
T2_EXTENSION_LIMIT_RC = 75  # Extension limit for reverse complement
T2_MISMATCH_THRESHOLD_EARLY = 10  # Early mismatch threshold
T2_MISMATCH_RATE = 0.08  # Mismatch rate threshold
T2_GAP_FILL_START_MAX = 60  # Maximum gap size for start extension
T2_GAP_FILL_BETWEEN_MAX = 40  # Maximum gap size for between-segment filling
T2_GAP_FILL_END_MAX = 60  # Maximum gap size for end extension

# T1 parameters (long sequences) - current best performing values
T1_EXTENSION_LIMIT_LONG = 198  # Extension limit for long sequences
T1_EXTENSION_LIMIT_RC = 198  # Extension limit for reverse complement
T1_MISMATCH_THRESHOLD_EARLY = 15  # Early mismatch threshold
T1_MISMATCH_RATE = 0.09  # Mismatch rate threshold
T1_GAP_FILL_START_MAX = 50  # Maximum gap size for start extension
T1_GAP_FILL_BETWEEN_MAX = 30  # Maximum gap size for between-segment filling
T1_GAP_FILL_END_MAX = 50  # Maximum gap size for end extension
