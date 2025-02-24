query = reference = ''

with open('Duplicate identification/input/query.txt') as f:
    query = f.read()

with open('Duplicate identification/input/reference.txt') as f:
    reference = f.read()

""""from genome_utils import inv

def find_duplicates(reference, query):
    results = []
    i = 0
    while i < len(query):
        match_size = 0
        # try to extend the match at position i
        for size in range(1, len(query) - i + 1):
            segment = query[i:i+size]
            if segment in reference or inv(segment) in reference:
                match_size = size
            else:
                break
        if match_size > 0:
            segment = query[i:i+match_size]
            is_inverse = segment not in reference
            # For simplicity, treat each match as count = 1
            results.append((i, match_size, 1, is_inverse))
            i += match_size
        else:
            i += 1
    return results

duplicates = find_duplicates(reference, query)
print(duplicates)

