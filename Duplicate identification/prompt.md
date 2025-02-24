# Find Duplicates of bp sequences in query sequence according to reference sequence

## Description

Find Duplicates of bp sequences in query sequence according to reference sequence.

We consider two sequences to be duplicates if one of the two conditions is met:
1. The two sequences are identical
2. The `inv()` of one of the sequences is identical to the other sequence

```python
BP_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def inv(s: str):
    return ''.join([BP_MAP[c] for c in s[::-1]])
```

### Example

s1, s2, s3 are all genome sequences.
reference sequence: s1 s2 s3
query sequence: s1 s2 s3 s3 inv(s2 s3) s3 

Then, s3 duplicates once, and s2 s3 duplicates once.

### 约定

The query sequence is formed from randomly inserting duplicates of the substring of the reference sequence into some position of reference sequence itself. 

If two substrings of reference sequence both have duplicates in query sequence and one substring involves the other within it, we only count the longer one.

### Goal

Implement an efficient algorithm for identifying duplicates within the provided sequences.

## Input & Output

Input: reference sequence, query sequence
Output: The starting index, size and consecutive duplicate count of each duplicate sequence in reference sequence, and if inversed duplicate of not.

genome_utils.py:
```python
BP_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def inv(s: str):
    return ''.join([BP_MAP[c] for c in s[::-1]])
```

This is main.py, but this program has some bugs.
```python
from genome_utils import inv

class State:
    def __init__(self):
        self.next = dict()
        self.link = -1
        self.len = 0

class SuffixAutomaton:
    def __init__(self):
        self.size = 1
        self.last = 0
        self.states = [State()]
    
    def sa_extend(self, c):
        p = self.last
        curr = self.size
        self.size += 1
        self.states.append(State())
        self.states[curr].len = self.states[p].len + 1
        while p >= 0 and c not in self.states[p].next:
            self.states[p].next[c] = curr
            p = self.states[p].link
        if p == -1:
            self.states[curr].link = 0
        else:
            q = self.states[p].next[c]
            if self.states[p].len + 1 == self.states[q].len:
                self.states[curr].link = q
            else:
                clone = self.size
                self.size += 1
                self.states.append(State())
                self.states[clone].len = self.states[p].len + 1
                self.states[clone].next = self.states[q].next.copy()
                self.states[clone].link = self.states[q].link
                while p >= 0 and self.states[p].next[c] == q:
                    self.states[p].next[c] = clone
                    p = self.states[p].link
                self.states[q].link = clone
                self.states[curr].link = clone
        self.last = curr

def build_suffix_automaton(s):
    sam = SuffixAutomaton()
    for c in s:
        sam.sa_extend(c)
    return sam

def find_max_length(sam, query, start):
    current_state = 0
    max_len = 0
    for j in range(start, len(query)):
        c = query[j]
        if c in sam.states[current_state].next:
            current_state = sam.states[current_state].next[c]
            max_len += 1
        else:
            break
    return max_len

# Read input
with open('Duplicate identification/input/reference.txt') as f:
    reference = f.read().strip('\n')

with open('Duplicate identification/input/query.txt') as f:
    query = f.read().strip('\n')

# Preprocess inv_reference
inv_reference = inv(reference)

# Build suffix automata for reference and inv_reference
sam_ref = build_suffix_automaton(reference)
sam_inv_ref = build_suffix_automaton(inv_reference)

# Precompute max lengths for each position in query
max_lengths = []
for i in range(len(query)):
    len_ref = find_max_length(sam_ref, query, i)
    len_inv = find_max_length(sam_inv_ref, query, i)
    max_len = max(len_ref, len_inv)
    inverted = len_inv > len_ref or (len_inv == len_ref and len_inv > 0)
    max_lengths.append( (max_len, inverted) )

# Process query to find all valid duplicates, considering maximum lengths and non-overlapping
result = []
i = 0
while i < len(query):
    max_len, is_inverted = max_lengths[i]
    if max_len == 0:
        i += 1
        continue
    # Determine the actual substring and its reference occurrence
    substr = query[i:i+max_len]
    if is_inverted:
        ref_substr = inv(substr)
    else:
        ref_substr = substr
    # Find the first occurrence in reference
    ref_start = reference.find(ref_substr)
    if ref_start == -1:
        i += max_len
        continue
    
    # Check consecutive duplicates
    consecutive = 1
    next_i = i + max_len
    while next_i <= len(query) - max_len:
        # Check if the next substring is the same as the current
        next_substr = query[next_i:next_i+max_len]
        if next_substr != substr:
            break
        # Check if the next position's max_len is the same and the inversion status matches
        next_max_len, next_is_inverted = max_lengths[next_i]
        if next_max_len != max_len or next_is_inverted != is_inverted:
            break
        # Check if the reference substring is the same
        current_ref_substr = inv(next_substr) if next_is_inverted else next_substr
        if current_ref_substr != ref_substr:
            break
        consecutive += 1
        next_i += max_len
    # Append result
    result.append( (ref_start, max_len, consecutive, is_inverted) )
    # Move to next position
    i = next_i

# Output the result
for entry in result:
    print(f"{entry[0]} {entry[1]} {entry[2]} {entry[3]}")
```

This program's output is:

0 402 1 False
352 50 3 False
352 48 1 False
330 70 3 False
298 102 1 True
300 98 1 True
400 400 1 False

while the expected output is:

350 50 4 False
330 70 3 False
300 100 2 True

Analyze what caused the bug and fix it.
