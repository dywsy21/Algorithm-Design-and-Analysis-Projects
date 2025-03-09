BP_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def inv(s: str):
    return ''.join([BP_MAP[c] for c in s[::-1]])


