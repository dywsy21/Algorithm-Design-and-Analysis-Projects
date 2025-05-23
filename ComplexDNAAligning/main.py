from evaluate import calculate_value, seq2hashtable_multi_test, calculate_distance
import config

def find_alignment_tuples(query, ref):
    # Generate anchors with multiple k-mer sizes for better coverage
    all_anchors = []
    for kmer_size in [12, 15, 18]:  # Multiple k-mer sizes
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, 1)  # Use shift=1 for maximum coverage
        all_anchors.extend(anchors)
    
    if len(all_anchors) == 0:
        return []
    
    # Extend anchors conservatively to ensure validation passes
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor_conservatively(query, ref, query_pos, ref_pos, strand, match_len)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            segments.append(segment)
    
    # Remove overlaps and select best coverage - FIXED
    return select_non_overlapping_segments(segments, query, ref)

def extend_anchor_conservatively(query, ref, qpos, rpos, strand, initial_len):
    """Conservatively extend anchor ensuring segments will pass validation"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:  # Forward strand
        # Extend left with exact matching only
        left_ext = 0
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and 
               query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]):
            left_ext += 1
            if left_ext > 500:  # Reasonable limit
                break
        
        # Extend right with exact matching only
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and
               query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]):
            right_ext += 1
            if right_ext > 500:  # Reasonable limit
                break
        
        qstart = qpos - left_ext
        qend = qpos + initial_len + right_ext
        rstart = rpos - left_ext
        rend = rpos + initial_len + right_ext
        
        return (qstart, qend, rstart, rend)
    
    else:  # Reverse complement strand
        # Extend left (query left, ref right) with exact matching
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
            if left_ext > 500:
                break
        
        # Extend right (query right, ref left) with exact matching
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
            if right_ext > 500:
                break
        
        qstart = qpos - left_ext
        qend = qpos + initial_len + right_ext
        rstart = rpos - right_ext
        rend = rpos + initial_len + left_ext
        
        return (qstart, qend, rstart, rend)

def select_non_overlapping_segments(segments, query, ref):
    """Select non-overlapping segments ensuring no query position overlap"""
    if not segments:
        return []
    
    # Remove duplicate segments
    segments = list(set(segments))
    
    # Filter out invalid segments
    valid_segments = []
    for seg in segments:
        qstart, qend, rstart, rend = seg
        if (qend - qstart >= config.MIN_SEGMENT_LENGTH and 
            qend <= len(query) and rend <= len(ref) and
            qstart >= 0 and rstart >= 0 and qend > qstart and rend > rstart):
            valid_segments.append(seg)
    
    if not valid_segments:
        return []
    
    # Sort by query start position
    valid_segments.sort(key=lambda x: x[0])
    
    # Use greedy algorithm to select non-overlapping segments with maximum total length
    selected = []
    last_query_end = 0
    
    for segment in valid_segments:
        qstart, qend, rstart, rend = segment
        
        # No overlap with previous segment
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    
    # Try to fill gaps between selected segments
    return fill_gaps_between_segments(selected, query, ref)

def fill_gaps_between_segments(segments, query, ref):
    """Fill gaps between segments where possible"""
    if len(segments) <= 1:
        return segments
    
    filled = []
    
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        
        if i == 0:
            # Try to extend first segment to start from 0
            if qstart > 0:
                gap_size = qstart
                if gap_size <= 100 and rstart >= gap_size:  # Conservative gap filling
                    # Check if we can extend with exact match
                    can_extend = True
                    for j in range(gap_size):
                        if query[j] != ref[rstart - gap_size + j]:
                            can_extend = False
                            break
                    
                    if can_extend:
                        qstart = 0
                        rstart = rstart - gap_size
        
        # Try to fill gap with previous segment
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            
            if 0 < gap_size <= 50:  # Small gap
                # Try to extend previous segment
                prev_segment = filled[-1]
                can_extend = True
                for j in range(gap_size):
                    if (prev_qend + j >= len(query) or 
                        prev_segment[3] + j >= len(ref) or
                        query[prev_qend + j] != ref[prev_segment[3] + j]):
                        can_extend = False
                        break
                
                if can_extend:
                    # Merge with previous segment
                    filled[-1] = (prev_segment[0], qend, prev_segment[2], rend)
                    continue
        
        filled.append((qstart, qend, rstart, rend))
    
    # Try to extend last segment to end of query
    if filled and filled[-1][1] < len(query):
        last_segment = filled[-1]
        gap_size = len(query) - last_segment[1]
        if gap_size <= 100 and last_segment[3] + gap_size <= len(ref):
            # Check if we can extend with exact match
            can_extend = True
            for j in range(gap_size):
                if query[last_segment[1] + j] != ref[last_segment[3] + j]:
                    can_extend = False
                    break
            
            if can_extend:
                filled[-1] = (last_segment[0], len(query), last_segment[2], last_segment[3] + gap_size)
    
    # Final validation: ensure no overlaps
    final_result = []
    last_end = 0
    for segment in filled:
        if segment[0] >= last_end:
            final_result.append(segment)
            last_end = segment[1]
    
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
