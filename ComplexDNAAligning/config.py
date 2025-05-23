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

# Extension parameters
INITIAL_EXTEND = 100     # Try to extend this much initially
MAX_MISMATCH_RATE = 0.05 # Allow small mismatch rate during extension
CHUNK_SIZE = 20          # Process in chunks to balance speed vs accuracy
MAX_EXTEND_LENGTH = 1000
MAX_GAP_FILL = 200

# Short sequence specific parameters
SHORT_SEQUENCE_THRESHOLD = 3000
SHORT_KMER_SIZES = [4, 6, 8, 10, 12]  # Multiple small k-mer sizes for sensitivity
SHORT_MAX_GAP_FILL = 300  # More aggressive gap filling for short sequences
SHORT_MISMATCH_TOLERANCE = 0.12  # Higher tolerance for short sequences
SHORT_EXTENSION_LIMIT = 150  # Conservative extension limit
