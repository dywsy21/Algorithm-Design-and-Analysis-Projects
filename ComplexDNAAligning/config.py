# Core parameters
KMER_SIZE = 15
MIN_SEGMENT_LENGTH = 30  # Required by scoring function
MAX_EDIT_RATE = 0.1      # Required by scoring function

# Multi-strategy parameters
SMALL_KMER_SIZE = 8      # For sensitive detection
LARGE_KMER_SIZE = 20     # For high-confidence anchors
DENSE_SHIFT = 1          # Dense sampling for better coverage
SPARSE_SHIFT = 3         # Sparse sampling for speed

# Sequence length threshold for strategy selection
SHORT_SEQUENCE_THRESHOLD = 3000

# T2 parameters (short sequences) - current best performing values
T2_KMER_SIZES = [6, 7, 8, 9, 11, 13, 15]  # K-mer sizes for short sequences
T2_EXTENSION_LIMIT_SHORT = 74  # Extension limit for short sequences
T2_EXTENSION_LIMIT_RC = 75  # Extension limit for reverse complement
T2_MISMATCH_THRESHOLD_EARLY = 11  # Early mismatch threshold
T2_MISMATCH_RATE = 0.08  # Mismatch rate threshold
T2_GAP_FILL_START_MAX = 60  # Maximum gap size for start extension
T2_GAP_FILL_BETWEEN_MAX = 40  # Maximum gap size for between-segment filling
T2_GAP_FILL_END_MAX = 60  # Maximum gap size for end extension

# T1 parameters (long sequences) - current best performing values
T1_EXTENSION_LIMIT_LONG = 200  # Extension limit for long sequences
T1_EXTENSION_LIMIT_RC = 200  # Extension limit for reverse complement
T1_MISMATCH_THRESHOLD_EARLY = 15  # Early mismatch threshold
T1_MISMATCH_RATE = 0.05  # Mismatch rate threshold
T1_GAP_FILL_START_MAX = 50  # Maximum gap size for start extension
T1_GAP_FILL_BETWEEN_MAX = 30  # Maximum gap size for between-segment filling
T1_GAP_FILL_END_MAX = 50  # Maximum gap size for end extension
