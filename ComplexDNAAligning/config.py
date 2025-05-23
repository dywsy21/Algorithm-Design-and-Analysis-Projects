# Core parameters
KMER_SIZE = 15
SHIFT = 5
MIN_SEGMENT_LENGTH = 30  # Required by scoring function
MAX_EDIT_RATE = 0.1      # Required by scoring function

# Extension parameters - use exact match extension first
INITIAL_EXTEND = 100     # Try to extend this much initially
MAX_MISMATCH_RATE = 0.05 # Allow small mismatch rate during extension
CHUNK_SIZE = 20          # Process in chunks to balance speed vs accuracy
