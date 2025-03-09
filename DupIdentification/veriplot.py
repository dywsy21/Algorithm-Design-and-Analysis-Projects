from genome_utils import *

# New function to verify duplicates and plot duplicate blocks (existing functionality)
def verify(result, query, reference):
    verified = True
    for (q_idx, ref_start, dup_len, count, inverted) in result:
        expected = query[q_idx:q_idx+dup_len*count]
        orig = reference[ref_start:ref_start+dup_len]
        if inverted:
            orig = inv(orig)
        if expected != orig * count:
            print(f"Verification fails at query index {q_idx}")
            verified = False
    if verified:
        print("All duplicates verified successfully.")

# New function to plot duplicate mapping:
def plot_duplicate_mapping(result, query, reference):
    import matplotlib.pyplot as plt
    mapping_x = []
    mapping_y = []
    # For each duplicate block, compute mapping for each base
    for (q_idx, ref_start, dup_len, count, inverted) in result:
        for occ in range(count):
            # For each occurrence in query:
            for b in range(dup_len):
                # Query index:
                q_pos = q_idx + occ * dup_len + b
                # For non-inverted, reference mapping is straightforward;
                # for inverted, the first base of query corresponds to last base of reference block.
                if not inverted:
                    r_pos = ref_start + b
                else:
                    r_pos = ref_start + (dup_len - 1 - b)
                mapping_x.append(q_pos)
                mapping_y.append(r_pos)
    plt.figure(figsize=(10,4))
    plt.scatter(mapping_x, mapping_y, s=10, color='tab:red')
    plt.xlabel("Query Sequence Index")
    plt.ylabel("Reference Sequence Index")
    plt.title("Mapping of Duplicate Bases from Query to Reference")
    plt.show()

