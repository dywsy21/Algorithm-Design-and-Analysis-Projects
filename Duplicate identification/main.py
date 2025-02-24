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
    if max_lengths[i][0] == 0:
        i += 1
        continue
    dup_len = max_lengths[i][0]
    current_inverted = max_lengths[i][1]
    start_sub = query[i:i+dup_len]
    count = 0
    j = i
    while j <= len(query) - dup_len:
        # Instead of shrinking the duplicate length, require every duplicate to have at least dup_len characters.
        if max_lengths[j][0] < dup_len or max_lengths[j][1] != current_inverted:
            break
        if query[j:j+dup_len] != start_sub:
            break
        count += 1
        j += dup_len
    ref_substr = inv(start_sub) if current_inverted else start_sub
    ref_start = reference.find(ref_substr)
    result.append((ref_start, dup_len, count, current_inverted))
    i = j

# Output the result
for entry in result:
    print(f"{entry[0]} {entry[1]} {entry[2]} {entry[3]}")
