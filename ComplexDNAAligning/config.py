# Core parameters
KMER_SIZE = 15
SHIFT = 5
MIN_SEGMENT_LENGTH = 30  # Required by scoring function
MAX_EDIT_RATE = 0.1      # Required by scoring function

# Multi-strategy parameters
SMALL_KMER_SIZE = 8      # For sensitive detection
LARGE_KMER_SIZE = 20     # For high-confidence anchors
DENSE_SHIFT = 1          # Dense sampling for better coverage
SPARSE_SHIFT = 3         # Sparse sampling for speed

# Sequence length threshold for strategy selection
SHORT_SEQUENCE_THRESHOLD = 3000
