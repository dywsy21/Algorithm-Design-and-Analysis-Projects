from evaluate import calculate_value, seq2hashtable_multi_test, calculate_distance
import config

def find_alignment_tuples(query, ref):
    """Main entry point - determine strategy based on sequence length"""
    if len(query) < config.SHORT_SEQUENCE_THRESHOLD:
        return find_alignment_with_strategy(query, ref, "short")
    else:
        return find_alignment_with_strategy(query, ref, "long")

def find_alignment_with_strategy(query, ref, strategy):
    """Unified alignment function with different strategies"""
    # Get strategy-specific parameters
    params = get_strategy_params(strategy)
    
    # Generate anchors
    all_anchors = []
    for kmer_size in params['kmer_sizes']:
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, params['shift'])
        all_anchors.extend(anchors)
    
    if not all_anchors:
        return []
    
    # Extend anchors
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor(query, ref, query_pos, ref_pos, strand, match_len, params)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            segments.append(segment)
    
    # Select and optimize segments
    return select_and_optimize_segments(segments, query, ref, params)

def get_strategy_params(strategy):
    """Get parameters for different strategies"""
    if strategy == "short":
        return {
            'kmer_sizes': [6, 8, 10, 12, 15],
            'shift': 1,
            'extension_limit': 100,
            'mismatch_tolerance': 0.08,
            'gap_fill_limit': 100,
            'merge_enabled': True,
            'aggressive_gap_fill': True
        }
    else:  # long
        return {
            'kmer_sizes': [config.SMALL_KMER_SIZE, config.KMER_SIZE],
            'shift': config.DENSE_SHIFT,
            'extension_limit': 200,
            'mismatch_tolerance': 0.05,
            'gap_fill_limit': 50,
            'merge_enabled': False,
            'aggressive_gap_fill': False
        }

def extend_anchor(query, ref, qpos, rpos, strand, initial_len, params):
    """Unified anchor extension with strategy-specific parameters"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:  # Forward strand
        # Left extension
        left_ext = 0
        mismatches = 0
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and 
               left_ext < params['extension_limit']):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
            else:
                mismatches += 1
                if (left_ext > 15 and 
                    mismatches / left_ext > params['mismatch_tolerance']):
                    break
                left_ext += 1
        
        # Right extension
        right_ext = 0
        mismatches = 0
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < params['extension_limit']):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
            else:
                mismatches += 1
                if (right_ext > 15 and 
                    mismatches / right_ext > params['mismatch_tolerance']):
                    break
                right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    
    else:  # Reverse complement strand - exact matching only
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               left_ext < params['extension_limit'] and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
        
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               right_ext < params['extension_limit'] and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def select_and_optimize_segments(segments, query, ref, params):
    """Unified segment selection and optimization"""
    # Validate segments
    valid_segments = []
    for seg in list(set(segments)):  # Remove duplicates
        qstart, qend, rstart, rend = seg
        if (qend - qstart >= config.MIN_SEGMENT_LENGTH and 
            qend <= len(query) and rend <= len(ref) and
            qstart >= 0 and rstart >= 0 and qend > qstart and
            is_segment_valid(query, ref, qstart, qend, rstart, rend)):
            valid_segments.append(seg)
    
    if not valid_segments:
        return []
    
    # Sort by query start position
    valid_segments.sort(key=lambda x: x[0])
    
    # Select non-overlapping segments
    selected = select_non_overlapping(valid_segments, params)
    
    # Apply gap filling based on strategy
    if params['aggressive_gap_fill']:
        return aggressive_gap_filling(selected, query, ref, params)
    else:
        return conservative_gap_filling(selected, query, ref, params)

def select_non_overlapping(segments, params):
    """Select non-overlapping segments with optional merging"""
    selected = []
    last_query_end = 0
    
    for segment in segments:
        qstart, qend, rstart, rend = segment
        
        # Try merging if enabled and there's a small gap
        if (params['merge_enabled'] and selected and 
            qstart <= last_query_end + 50):
            merged = try_merge_segments(selected[-1], segment)
            if merged:
                selected[-1] = merged
                last_query_end = merged[1]
                continue
        
        # Add if no overlap
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    
    return selected

def try_merge_segments(seg1, seg2):
    """Try to merge two segments with basic validation"""
    q1_start, q1_end, r1_start, r1_end = seg1
    q2_start, q2_end, r2_start, r2_end = seg2
    
    q_gap = q2_start - q1_end
    r_gap = abs(r2_start - r1_end)
    
    if q_gap < 0 or q_gap > 100 or r_gap > 200:
        return None
    
    # Simple merge - just span the coordinates
    merged_q_start = q1_start
    merged_q_end = q2_end
    merged_r_start = min(r1_start, r2_start)
    merged_r_end = max(r1_end, r2_end)
    
    return (merged_q_start, merged_q_end, merged_r_start, merged_r_end)

def aggressive_gap_filling(segments, query, ref, params):
    """Aggressive gap filling for short sequences"""
    if not segments:
        return []
    
    filled = list(segments)
    qlen = len(query)
    
    # Extend first segment to start
    if filled[0][0] > 0:
        first_seg = filled[0]
        gap_size = first_seg[0]
        if gap_size <= params['gap_fill_limit'] and first_seg[2] >= gap_size:
            extended = (0, first_seg[1], first_seg[2] - gap_size, first_seg[3])
            if is_segment_valid(query, ref, *extended):
                filled[0] = extended
    
    # Extend last segment to end
    if filled[-1][1] < qlen:
        last_seg = filled[-1]
        gap_size = qlen - last_seg[1]
        if gap_size <= params['gap_fill_limit'] and last_seg[3] + gap_size <= len(ref):
            extended = (last_seg[0], qlen, last_seg[2], last_seg[3] + gap_size)
            if is_segment_valid(query, ref, *extended):
                filled[-1] = extended
    
    # Try to merge adjacent segments
    merged = []
    for segment in filled:
        if merged:
            merged_seg = try_merge_segments(merged[-1], segment)
            if merged_seg and is_segment_valid(query, ref, *merged_seg):
                merged[-1] = merged_seg
                continue
        merged.append(segment)
    
    return validate_no_overlaps(merged)

def conservative_gap_filling(segments, query, ref, params):
    """Conservative gap filling for long sequences"""
    if not segments:
        return []
    
    filled = []
    qlen = len(query)
    
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        
        # Extend first segment to start (small gaps only)
        if i == 0 and qstart > 0:
            gap_size = qstart
            if gap_size <= params['gap_fill_limit'] and rstart >= gap_size:
                # Check for exact match
                can_extend = all(query[j] == ref[rstart - gap_size + j] 
                               for j in range(gap_size))
                if can_extend:
                    qstart, rstart = 0, rstart - gap_size
        
        # Fill small gaps between segments
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            
            if 0 < gap_size <= 30:  # Very small gaps
                prev_segment = filled[-1]
                can_extend = all(
                    prev_qend + j < len(query) and 
                    prev_segment[3] + j < len(ref) and
                    query[prev_qend + j] == ref[prev_segment[3] + j]
                    for j in range(gap_size)
                )
                
                if can_extend:
                    filled[-1] = (prev_segment[0], qend, prev_segment[2], rend)
                    continue
        
        filled.append((qstart, qend, rstart, rend))
    
    # Extend last segment to end
    if filled and filled[-1][1] < qlen:
        last_segment = filled[-1]
        gap_size = qlen - last_segment[1]
        if gap_size <= params['gap_fill_limit'] and last_segment[3] + gap_size <= len(ref):
            can_extend = all(query[last_segment[1] + j] == ref[last_segment[3] + j]
                           for j in range(gap_size))
            if can_extend:
                filled[-1] = (last_segment[0], qlen, last_segment[2], last_segment[3] + gap_size)
    
    return validate_no_overlaps(filled)

def is_segment_valid(query, ref, qstart, qend, rstart, rend):
    """Validate segment using edit distance"""
    if qend - qstart < config.MIN_SEGMENT_LENGTH:
        return False
    
    try:
        edit_dist = calculate_distance(ref, query, rstart, rend, qstart, qend)
        edit_rate = edit_dist / (qend - qstart)
        return edit_rate <= config.MAX_EDIT_RATE
    except:
        return False

def validate_no_overlaps(segments):
    """Final validation to ensure no overlapping query positions"""
    if not segments:
        return []
    
    segments.sort(key=lambda x: x[0])
    validated = []
    last_end = 0
    
    for segment in segments:
        qstart, qend, rstart, rend = segment
        if qstart >= last_end and qend - qstart >= config.MIN_SEGMENT_LENGTH:
            validated.append((qstart, qend, rstart, rend))
            last_end = qend
    
    return validated

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
