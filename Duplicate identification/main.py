query = reference = ''

with open('Duplicate identification/input/query.txt') as f:
    query = f.read().strip('\n')
    assert all(c in 'ATCG' for c in query)

with open('Duplicate identification/input/reference.txt') as f:
    reference = f.read().strip('\n')
    assert all(c in 'ATCG' for c in reference)

from genome_utils import inv

def find_duplicates(reference, query):
    results = []
    i = 0
    while i < len(query):
        match_size = 0
        is_inverse = False
        ref_index = -1
        
        for size in range(1, len(query) - i + 1):
            s = query[i:i+size]
            if s in reference:
                match_size = size
                is_inverse = False
                ref_index = reference.index(s)
            elif inv(s) in reference:
                match_size = size
                is_inverse = True
                ref_index = reference.index(inv(s))
            else:
                break
                
        if match_size > 0:
            count = 1
            j = i + match_size
            while j + match_size <= len(query):
                next_seg = query[j:j+match_size]
                if next_seg == query[i:i+match_size] or next_seg == inv(query[i:i+match_size]):
                    count += 1
                    j += match_size
                else:
                    break
            results.append((i, ref_index, match_size, count, is_inverse))
            i = j
        else:
            i += 1
    
    # prune repeating count == 1 results
    return [r for r in results if r[3] > 1]

def reconstruct_query(duplicates, reference):
    reconstructed = ''
    last_query_idx = 0
    
    # Sort by query index to reconstruct in order
    sorted_dups = sorted(duplicates, key=lambda x: x[0])
    
    for query_idx, ref_idx, size, count, is_inverse in sorted_dups:
        # Fill any gap before this match
        if query_idx > last_query_idx:
            missing_part = query[last_query_idx:query_idx]
            reconstructed += missing_part
            
        # Get the matching segment from reference
        segment = reference[ref_idx:ref_idx + size]
        if is_inverse:
            segment = inv(segment)
            
        # Add it count times
        reconstructed += segment * count
        last_query_idx = query_idx + (size * count)
    
    # Add any remaining part after the last match
    if last_query_idx < len(query):
        reconstructed += query[last_query_idx:]
        
    return reconstructed

print(reference.__len__(), query.__len__())
duplicates = find_duplicates(reference, query)
print("Duplicates found:", duplicates)

reconstructed = reconstruct_query(duplicates, reference)
print("Reconstruction successful:", reconstructed == query)
if reconstructed != query:
    print("Original length:", len(query))
    print("Reconstructed length:", len(reconstructed))
    # Find first difference
    for i, (a, b) in enumerate(zip(query, reconstructed)):
        if a != b:
            print(f"First difference at index {i}: {a} vs {b}")
            break
