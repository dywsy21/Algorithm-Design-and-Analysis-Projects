# Find Duplicates of bp sequences in query sequence according to reference sequence

## io

Input: reference sequence, query sequence
Output: The starting index, size and consecutive duplicate count of each duplicate sequence in reference sequence, and if inversed duplicate of not.

## Description

Find Duplicates of bp sequences in query sequence according to reference sequence.

We consider two sequences to be duplicates if one of the two conditions is met:
1. The two sequences are identical
2. The `get_inverse` of one of the sequences is identical to the other sequence

```python
BP_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def inv(s: str):
    return ''.join([BP_MAP[c] for c in s[::-1]])
```

For example: s1, s2, s3 are all genome sequences.
reference sequence: s1 s2 s3
query sequence: s1 s2 s3 s3 inv(s2 s3) s3 

Then, s3 duplicates once, and s2 s3 duplicates once.
