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
        for size in range(1, len(query) - i + 1):
            s = query[i:i+size]
            if s in reference:
                match_size = size
                is_inverse = False
            elif inv(s) in reference:
                match_size = size
                is_inverse = True
            else:
                break
        if match_size > 0:
            segment = query[i:i+match_size]
            count = 1
            j = i + match_size
            while j + match_size <= len(query):
                next_seg = query[j:j+match_size]
                if next_seg == segment or next_seg == inv(segment):
                    count += 1
                    j += match_size
                else:
                    break
            results.append((i, match_size, count, is_inverse))
            i = j
        else:
            i += 1
    return results

print(reference.__len__(), query.__len__())
duplicates = find_duplicates(reference, query)
print(duplicates)

