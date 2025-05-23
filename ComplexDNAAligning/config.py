# K-mer parameters
KMER_SIZE = 12  # Smaller k-mer for better sensitivity
SHIFT = 3       # Smaller shift for more anchors

# Extension parameters
MIN_ANCHOR_LENGTH = 25    # Minimum length for valid anchors
EXTEND_MISMATCH_RATE = 0.15  # Allow some mismatches during extension

# Chaining parameters
MAX_GAP_RATIO = 3.0      # Maximum gap ratio between anchors
MIN_CHAIN_SCORE = 50     # Minimum score to consider a chain
CHAIN_GAP_PENALTY = 0.1  # Penalty for gaps in chaining

# Gap filling parameters
FILL_SMALL_GAPS = True   # Whether to fill small gaps between anchors
MAX_FILL_GAP = 500       # Maximum gap size to attempt filling
