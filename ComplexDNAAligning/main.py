from evaluate import calculate_value, seq2hashtable_multi_test, calculate_distance
import config

def find_alignment_tuples(query, ref):
    # Determine if this is a short sequence and use appropriate strategy
    if len(query) < config.SHORT_SEQUENCE_THRESHOLD:
        return find_alignment_tuples_short_optimized(query, ref)
    else:
        return find_alignment_tuples_long(query, ref)

def find_alignment_tuples_short_optimized(query, ref):
    """Conservative but comprehensive strategy for short sequences"""
    # Use multiple k-mer sizes but be more conservative
    all_anchors = []
    
    # Focus on proven k-mer sizes that work well
    for kmer_size in [6, 8, 10, 12, 15]:
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, 1)
        all_anchors.extend(anchors)
    
    if len(all_anchors) == 0:
        return []
    
    # Conservative extension to ensure segments pass validation
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor_conservatively_short(query, ref, query_pos, ref_pos, strand, match_len)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            segments.append(segment)
    
    # Conservative selection with intelligent merging
    return select_and_merge_short_sequences(segments, query, ref)

def find_alignment_tuples_long(query, ref):
    """Proven strategy for long sequences"""
    # Generate anchors with different k-mer sizes
    all_anchors = []
    
    # Use fewer strategies to reduce overlap complexity
    for kmer_size in [config.SMALL_KMER_SIZE, config.KMER_SIZE]:
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, config.DENSE_SHIFT)
        all_anchors.extend(anchors)
    
    if len(all_anchors) == 0:
        return []
    
    # Extend anchors conservatively
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor_conservatively(query, ref, query_pos, ref_pos, strand, match_len)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            segments.append(segment)
    
    # Ensure no overlapping segments
    return select_non_overlapping_greedy(segments, query, ref)

def extend_anchor_aggressively_short(query, ref, qpos, rpos, strand, initial_len):
    """Aggressive extension specifically optimized for short sequences"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:  # Forward strand
        # Aggressive left extension with higher mismatch tolerance
        left_ext = 0
        mismatches = 0
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and left_ext < config.SHORT_EXTENSION_LIMIT):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
            else:
                mismatches += 1
                # Higher tolerance for short sequences
                if left_ext > 10 and mismatches / left_ext > config.SHORT_MISMATCH_TOLERANCE:
                    break
                left_ext += 1
        
        # Aggressive right extension
        right_ext = 0
        mismatches = 0
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < config.SHORT_EXTENSION_LIMIT):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
            else:
                mismatches += 1
                if right_ext > 10 and mismatches / right_ext > config.SHORT_MISMATCH_TOLERANCE:
                    break
                right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    
    else:  # Reverse complement strand
        # Conservative reverse complement extension
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               left_ext < config.SHORT_EXTENSION_LIMIT and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
        
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               right_ext < config.SHORT_EXTENSION_LIMIT and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def extend_anchor_conservatively_short(query, ref, qpos, rpos, strand, initial_len):
    """Conservative extension for short sequences ensuring validation passes"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:  # Forward strand
        # Conservative left extension - prioritize exact matches
        left_ext = 0
        mismatches = 0
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and left_ext < 100):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
            else:
                mismatches += 1
                # Very strict for validation passing
                if left_ext > 15 and mismatches / left_ext > 0.08:
                    break
                left_ext += 1
        
        # Conservative right extension
        right_ext = 0
        mismatches = 0
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < 100):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
            else:
                mismatches += 1
                if right_ext > 15 and mismatches / right_ext > 0.08:
                    break
                right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    
    else:  # Reverse complement strand - exact matching only for reliability
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               left_ext < 100 and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
        
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               right_ext < 100 and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def select_and_merge_short_sequences(segments, query, ref):
    """Intelligent selection and merging for short sequences"""
    if not segments:
        return []
    
    # Remove duplicates and validate
    segments = list(set(segments))
    valid_segments = []
    for seg in segments:
        qstart, qend, rstart, rend = seg
        if (qend - qstart >= config.MIN_SEGMENT_LENGTH and 
            qend <= len(query) and rend <= len(ref) and
            qstart >= 0 and rstart >= 0 and qend > qstart):
            # Pre-validate with actual edit distance check
            if is_segment_valid(query, ref, qstart, qend, rstart, rend):
                valid_segments.append(seg)
    
    if not valid_segments:
        return []
    
    # Sort by query start position
    valid_segments.sort(key=lambda x: x[0])
    
    # Greedy selection with intelligent merging opportunities
    selected = []
    last_query_end = 0
    
    for segment in valid_segments:
        qstart, qend, rstart, rend = segment
        
        # Check for merging opportunity with the last segment
        if selected and qstart <= last_query_end + 50:  # Small gap or overlap
            last_segment = selected[-1]
            merged = try_merge_segments(last_segment, segment, query, ref)
            if merged:
                selected[-1] = merged
                last_query_end = merged[1]
                continue
        
        # No overlap with previous segment
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    
    # Final gap filling for short sequences
    return intelligent_gap_filling_short(selected, query, ref)

def is_segment_valid(query, ref, qstart, qend, rstart, rend):
    """Pre-validate segment using actual edit distance calculation"""
    if qend - qstart < config.MIN_SEGMENT_LENGTH:
        return False
    
    try:
        edit_dist = calculate_distance(ref, query, rstart, rend, qstart, qend)
        edit_rate = edit_dist / (qend - qstart)
        return edit_rate <= config.MAX_EDIT_RATE
    except:
        return False

def try_merge_segments(seg1, seg2, query, ref):
    """Try to merge two segments if the result is valid"""
    q1_start, q1_end, r1_start, r1_end = seg1
    q2_start, q2_end, r2_start, r2_end = seg2
    
    # Check if segments can be reasonably merged
    q_gap = q2_start - q1_end
    r_gap = abs(r2_start - r1_end)
    
    # Only try merging if gaps are reasonable
    if q_gap < 0 or q_gap > 100 or r_gap > 200:
        return None
    
    # Try to create merged segment
    merged_q_start = q1_start
    merged_q_end = q2_end
    
    # Determine reference coordinates based on strand direction
    if r1_end <= r2_start:  # Forward direction
        merged_r_start = r1_start
        merged_r_end = r2_end
    elif r2_end <= r1_start:  # Reverse direction
        merged_r_start = r2_start
        merged_r_end = r1_end
    else:
        # Overlapping ref coordinates - use the span
        merged_r_start = min(r1_start, r2_start)
        merged_r_end = max(r1_end, r2_end)
    
    # Validate the merged segment
    if is_segment_valid(query, ref, merged_q_start, merged_q_end, merged_r_start, merged_r_end):
        return (merged_q_start, merged_q_end, merged_r_start, merged_r_end)
    
    return None

def intelligent_gap_filling_short(segments, query, ref):
    """Intelligent gap filling for short sequences"""
    if not segments:
        return []
    
    filled = []
    qlen = len(query)
    
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        
        # Try to extend first segment to start from 0
        if i == 0 and qstart > 0:
            gap_size = qstart
            if gap_size <= 100 and rstart >= gap_size:
                extended_segment = (0, qend, rstart - gap_size, rend)
                if is_segment_valid(query, ref, *extended_segment):
                    qstart, rstart = 0, rstart - gap_size
        
        # Try to fill gaps between segments
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            
            if 0 < gap_size <= 80:  # Reasonable gap size
                prev_segment = filled[-1]
                merged_segment = (prev_segment[0], qend, prev_segment[2], rend)
                if is_segment_valid(query, ref, *merged_segment):
                    # Merge with previous segment
                    filled[-1] = merged_segment
                    continue
        
        filled.append((qstart, qend, rstart, rend))
    
    # Try to extend last segment to end of query
    if filled and filled[-1][1] < qlen:
        last_segment = filled[-1]
        gap_size = qlen - last_segment[1]
        if gap_size <= 100 and last_segment[3] + gap_size <= len(ref):
            extended_segment = (last_segment[0], qlen, last_segment[2], last_segment[3] + gap_size)
            if is_segment_valid(query, ref, *extended_segment):
                filled[-1] = extended_segment
    
    return filled

def extend_anchor_conservatively(query, ref, qpos, rpos, strand, initial_len):
    """Conservatively extend anchor with limited mismatch tolerance"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:  # Forward strand
        # Extend left conservatively
        left_ext = 0
        mismatches = 0
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and left_ext < 200):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
            else:
                mismatches += 1
                if left_ext > 10 and mismatches / left_ext > 0.05:  # Very conservative
                    break
                left_ext += 1
        
        # Extend right conservatively
        right_ext = 0
        mismatches = 0
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < 200):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
            else:
                mismatches += 1
                if right_ext > 10 and mismatches / right_ext > 0.05:
                    break
                right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    
    else:  # Reverse complement strand - exact matching only
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               left_ext < 200 and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
        
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               right_ext < 200 and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def select_non_overlapping_greedy(segments, query, ref):
    """Greedy selection ensuring absolutely no overlaps"""
    if not segments:
        return []
    
    # Remove duplicates and invalid segments
    segments = list(set(segments))
    valid_segments = []
    for seg in segments:
        qstart, qend, rstart, rend = seg
        if (qend - qstart >= config.MIN_SEGMENT_LENGTH and 
            qend <= len(query) and rend <= len(ref) and
            qstart >= 0 and rstart >= 0 and qend > qstart):
            valid_segments.append(seg)
    
    if not valid_segments:
        return []
    
    # Sort by query start position
    valid_segments.sort(key=lambda x: x[0])
    
    # Greedy selection with strict non-overlap enforcement
    selected = []
    last_query_end = 0
    
    for segment in valid_segments:
        qstart, qend, rstart, rend = segment
        
        # Strict non-overlap check
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    
    # Conservative gap filling
    return fill_gaps_conservatively(selected, query, ref)

def fill_gaps_conservatively(segments, query, ref):
    """Conservative gap filling to avoid creating overlaps"""
    if len(segments) <= 1:
        return segments
    
    filled = []
    qlen = len(query)
    
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        
        # Try to extend first segment to start from 0
        if i == 0 and qstart > 0:
            gap_size = qstart
            if gap_size <= 50 and rstart >= gap_size:  # Small gaps only
                # Check for exact match extension
                can_extend = True
                for j in range(gap_size):
                    if query[j] != ref[rstart - gap_size + j]:
                        can_extend = False
                        break
                
                if can_extend:
                    qstart = 0
                    rstart = rstart - gap_size
        
        # Try to fill small gaps between segments
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            
            if 0 < gap_size <= 30:  # Very small gaps only
                prev_segment = filled[-1]
                # Check for exact match extension
                can_extend = True
                for j in range(gap_size):
                    if (prev_qend + j >= len(query) or 
                        prev_segment[3] + j >= len(ref) or
                        query[prev_qend + j] != ref[prev_segment[3] + j]):
                        can_extend = False
                        break
                
                if can_extend:
                    # Merge segments
                    filled[-1] = (prev_segment[0], qend, prev_segment[2], rend)
                    continue
        
        filled.append((qstart, qend, rstart, rend))
    
    # Try to extend last segment to end of query
    if filled and filled[-1][1] < qlen:
        last_segment = filled[-1]
        gap_size = qlen - last_segment[1]
        if gap_size <= 50 and last_segment[3] + gap_size <= len(ref):
            # Check for exact match extension
            can_extend = True
            for j in range(gap_size):
                if query[last_segment[1] + j] != ref[last_segment[3] + j]:
                    can_extend = False
                    break
            
            if can_extend:
                filled[-1] = (last_segment[0], qlen, last_segment[2], last_segment[3] + gap_size)
    
    # Final validation: absolutely ensure no overlaps
    final_result = []
    last_end = 0
    for segment in filled:
        qstart, qend, rstart, rend = segment
        if qstart >= last_end and qend - qstart >= config.MIN_SEGMENT_LENGTH:
            final_result.append((qstart, qend, rstart, rend))
            last_end = qend
    
    return final_result

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
