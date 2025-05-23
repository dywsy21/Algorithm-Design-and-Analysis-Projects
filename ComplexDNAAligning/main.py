from evaluate import calculate_value, seq2hashtable_multi_test
import numpy as np
import config

def find_alignment_tuples(query, ref):
    # Step 1: Generate anchors with multiple k-mer sizes for better coverage
    all_anchors = []
    for kmer_size in [config.KMER_SIZE, config.KMER_SIZE + 3]:
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, config.SHIFT)
        if len(anchors) > 0:
            all_anchors.extend(anchors)
    
    if len(all_anchors) == 0:
        return []
    
    # Step 2: Extend anchors with mismatch tolerance
    extended_anchors = []
    for anchor in all_anchors:
        query_pos, ref_pos, strand, match_len = anchor
        extended_match = extend_match_with_mismatches(query, ref, query_pos, ref_pos, strand, match_len)
        if extended_match and extended_match[1] - extended_match[0] >= config.MIN_ANCHOR_LENGTH:
            extended_anchors.append(extended_match)
    
    if not extended_anchors:
        return []
    
    # Step 3: Remove redundant/overlapping anchors
    filtered_anchors = filter_overlapping_anchors(extended_anchors)
    
    # Step 4: Find optimal chain with gap penalties
    chain = find_optimal_chain_improved(filtered_anchors, query, ref)
    
    # Step 5: Fill gaps between anchors
    if config.FILL_SMALL_GAPS:
        chain = fill_gaps_between_anchors(chain, query, ref)
    
    return chain

def extend_match_with_mismatches(query, ref, query_pos, ref_pos, strand, initial_len):
    """Extend match allowing some mismatches"""
    query_len, ref_len = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:  # Forward
        # Extend left
        left_ext = 0
        mismatches = 0
        while (query_pos - left_ext > 0 and ref_pos - left_ext > 0):
            if query[query_pos - left_ext - 1] != ref[ref_pos - left_ext - 1]:
                mismatches += 1
                if mismatches / (left_ext + 1) > config.EXTEND_MISMATCH_RATE:
                    break
            left_ext += 1
            if left_ext > 100:  # Limit extension
                break
        
        # Extend right
        right_ext = 0
        mismatches = 0
        while (query_pos + initial_len + right_ext < query_len and 
               ref_pos + initial_len + right_ext < ref_len):
            if query[query_pos + initial_len + right_ext] != ref[ref_pos + initial_len + right_ext]:
                mismatches += 1
                if mismatches / (right_ext + 1) > config.EXTEND_MISMATCH_RATE:
                    break
            right_ext += 1
            if right_ext > 100:
                break
        
        return (query_pos - left_ext, query_pos + initial_len + right_ext,
                ref_pos - left_ext, ref_pos + initial_len + right_ext)
    else:  # Reverse complement
        # Similar logic for reverse complement
        left_ext = 0
        while (query_pos - left_ext > 0 and ref_pos + initial_len + left_ext < ref_len and
               query[query_pos - left_ext - 1] == rc_map.get(ref[ref_pos + initial_len + left_ext], 'N')):
            left_ext += 1
            if left_ext > 100:
                break
        
        right_ext = 0
        while (query_pos + initial_len + right_ext < query_len and ref_pos - right_ext > 0 and
               query[query_pos + initial_len + right_ext] == rc_map.get(ref[ref_pos - right_ext - 1], 'N')):
            right_ext += 1
            if right_ext > 100:
                break
        
        return (query_pos - left_ext, query_pos + initial_len + right_ext,
                ref_pos - right_ext, ref_pos + initial_len + left_ext)

def filter_overlapping_anchors(anchors):
    """Remove overlapping anchors, keeping the longer ones"""
    anchors = sorted(anchors, key=lambda x: x[1] - x[0], reverse=True)  # Sort by length
    filtered = []
    
    for anchor in anchors:
        overlap = False
        for existing in filtered:
            if (max(anchor[0], existing[0]) < min(anchor[1], existing[1]) and
                max(anchor[2], existing[2]) < min(anchor[3], existing[3])):
                overlap = True
                break
        if not overlap:
            filtered.append(anchor)
    
    return filtered

def find_optimal_chain_improved(anchors, query, ref):
    """Improved chaining with better scoring"""
    if not anchors:
        return []
    
    anchors = sorted(anchors, key=lambda x: x[0])
    n = len(anchors)
    dp = [0] * n
    prev = [-1] * n
    
    # Initialize with anchor lengths
    for i in range(n):
        dp[i] = anchors[i][1] - anchors[i][0]
    
    # DP with gap penalties
    for i in range(1, n):
        for j in range(i):
            if can_chain_improved(anchors[j], anchors[i]):
                gap_penalty = calculate_gap_penalty(anchors[j], anchors[i])
                score = dp[j] + (anchors[i][1] - anchors[i][0]) - gap_penalty
                if score > dp[i]:
                    dp[i] = score
                    prev[i] = j
    
    # Find best chain
    best_score = max(dp)
    if best_score < config.MIN_CHAIN_SCORE:
        return []
    
    best_end = dp.index(best_score)
    
    # Reconstruct chain
    chain = []
    current = best_end
    while current != -1:
        chain.append(anchors[current])
        current = prev[current]
    
    chain.reverse()
    return chain

def can_chain_improved(anchor1, anchor2):
    """Improved chaining logic"""
    q1_start, q1_end, r1_start, r1_end = anchor1
    q2_start, q2_end, r2_start, r2_end = anchor2
    
    if q2_start <= q1_end:
        return False
    
    gap_query = q2_start - q1_end
    gap_ref = abs(r2_start - r1_end)
    max_gap = max(gap_query * config.MAX_GAP_RATIO, config.MAX_FILL_GAP)
    
    return gap_ref <= max_gap

def calculate_gap_penalty(anchor1, anchor2):
    """Calculate penalty for gaps between anchors"""
    gap_query = anchor2[0] - anchor1[1]
    gap_ref = abs(anchor2[2] - anchor1[3])
    return max(gap_query, gap_ref) * config.CHAIN_GAP_PENALTY

def fill_gaps_between_anchors(chain, query, ref):
    """Fill small gaps between consecutive anchors"""
    if len(chain) <= 1:
        return chain
    
    filled_chain = [chain[0]]
    
    for i in range(1, len(chain)):
        prev_anchor = filled_chain[-1]
        curr_anchor = chain[i]
        
        gap_query = curr_anchor[0] - prev_anchor[1]
        gap_ref = curr_anchor[2] - prev_anchor[3]
        
        # Try to fill small gaps
        if 0 < gap_query <= config.MAX_FILL_GAP and 0 < gap_ref <= config.MAX_FILL_GAP:
            # Simple gap filling - just extend the previous anchor
            extended_anchor = (prev_anchor[0], curr_anchor[1], prev_anchor[2], curr_anchor[3])
            filled_chain[-1] = extended_anchor
        else:
            filled_chain.append(curr_anchor)
    
    return filled_chain

def main():
    q1 = q2 = r1 = r2 = ''
    with open('data/query1.txt', 'r') as f:
        q1 = f.read()
    with open('data/query2.txt', 'r') as f:
        q2 = f.read()
    with open('data/ref1.txt', 'r') as f:
        r1 = f.read()
    with open('data/ref2.txt', 'r') as f:
        r2 = f.read()
    q1 = q1.strip()
    q2 = q2.strip()
    r1 = r1.strip()
    r2 = r2.strip()
    t1 = find_alignment_tuples(q1, r1)
    t2 = find_alignment_tuples(q2, r2)
    print('t1: ', t1)
    print('t2: ', t2)
    print('t1 score: ', calculate_value(str(t1), r1, q1))
    print('t2 score: ', calculate_value(str(t2), r2, q2))


if __name__ == "__main__":
    main()
