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
        
        # Find the longest matching segment
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
            segment = query[i:i+match_size]
            count = 1
            j = i + match_size
            # Validate consecutive matches
            while j + match_size <= len(query):
                next_seg = query[j:j+match_size]
                current = segment if not is_inverse else inv(segment)
                if next_seg == current:
                    count += 1
                    j += match_size
                else:
                    break
            
            if count > 1:  # Only store duplicates
                results.append((i, ref_index, match_size, count, is_inverse))
            else:  # Store non-duplicate as regular segment
                results.append((i, ref_index, match_size, 1, is_inverse))
            i = j
        else:
            # Store unmatched single base as is
            results.append((i, -1, 1, 1, False))
            i += 1
    return results

def reconstruct_query(duplicates, reference):
    reconstructed = ''
    last_idx = 0
    
    for query_idx, ref_idx, size, count, is_inverse in sorted(duplicates, key=lambda x: x[0]):
        assert query_idx == last_idx, f"Gap in reconstruction at {last_idx} to {query_idx}"
        
        if ref_idx == -1:  # Unmatched base
            reconstructed += query[query_idx]  # This is the only place we need original query
        else:
            segment = reference[ref_idx:ref_idx + size]
            if is_inverse:
                segment = inv(segment)
            reconstructed += segment * count
            
        last_idx = query_idx + (size * count)
    
    return reconstructed

def main():
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

if __name__ == '__main__':
    main()
