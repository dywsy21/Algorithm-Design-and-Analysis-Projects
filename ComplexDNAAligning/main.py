from evaluate import calculate_value, seq2hashtable_multi_test, calculate_distance
import config

def find_alignment_tuples(query, ref):
    """Main entry point - determine strategy based on sequence length"""
    if len(query) < config.SHORT_SEQUENCE_THRESHOLD:
        return find_alignment_short_sequence(query, ref)
    else:
        return find_alignment_long_sequence(query, ref)

def find_alignment_short_sequence(query, ref):
    """Specialized strategy for short sequences (t2) - target 2090"""
    # Multi-pass anchor generation for better coverage
    all_anchors = []
    
    # First pass: standard k-mers
    for kmer_size in config.T2_KMER_SIZES:
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, 1)
        all_anchors.extend(anchors)
    
    # Second pass: fill potential gaps with very small k-mers
    if len(all_anchors) > 0:
        covered_regions = get_covered_regions(all_anchors, len(query))
        gap_regions = find_gap_regions(covered_regions, len(query))
        
        # Use smaller k-mers for gaps
        for gap_start, gap_end in gap_regions:
            if gap_end - gap_start > 50:  # Only for significant gaps
                gap_query = query[gap_start:gap_end]
                gap_anchors = seq2hashtable_multi_test(ref, gap_query, 4, 1)  # Very small k-mer
                # Adjust coordinates
                for qpos, rpos, strand, match_len in gap_anchors:
                    all_anchors.append((gap_start + qpos, rpos, strand, match_len))
    
    if not all_anchors:
        return []
    
    # Enhanced extension
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor_conservative_short(query, ref, query_pos, ref_pos, strand, match_len)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            if is_segment_valid(query, ref, *segment):
                segments.append(segment)
    
    return select_segments_conservative_short(segments, query, ref)

def find_alignment_long_sequence(query, ref):
    """Optimized strategy for long sequences (t1) - target 29820"""
    all_anchors = []
    
    # Standard anchor generation
    for kmer_size in [config.SMALL_KMER_SIZE, config.KMER_SIZE]:
        anchors = seq2hashtable_multi_test(ref, query, kmer_size, config.DENSE_SHIFT)
        all_anchors.extend(anchors)
    
    # Additional pass with larger k-mer for high specificity regions
    large_anchors = seq2hashtable_multi_test(ref, query, 20, 2)
    all_anchors.extend(large_anchors)
    
    if not all_anchors:
        return []
    
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor_conservative_long(query, ref, query_pos, ref_pos, strand, match_len)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            if is_segment_valid(query, ref, *segment):
                segments.append(segment)
    
    return select_segments_long(segments, query, ref)

def get_covered_regions(anchors, query_length):
    """Get regions covered by anchors"""
    covered = [False] * query_length
    for qpos, rpos, strand, match_len in anchors:
        for i in range(max(0, qpos - 10), min(query_length, qpos + match_len + 10)):
            covered[i] = True
    return covered

def find_gap_regions(covered, query_length):
    """Find uncovered gap regions"""
    gaps = []
    gap_start = None
    for i, is_covered in enumerate(covered):
        if not is_covered:
            if gap_start is None:
                gap_start = i
        else:
            if gap_start is not None:
                gaps.append((gap_start, i))
                gap_start = None
    
    if gap_start is not None:
        gaps.append((gap_start, query_length))
    
    return gaps

def extend_anchor_conservative_short(query, ref, qpos, rpos, strand, initial_len):
    """Conservative extension for short sequences with micro-optimizations"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:
        # Enhanced left extension with early stopping
        left_ext = 0
        mismatches = 0
        consecutive_mismatches = 0
        
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and 
               left_ext < config.T2_EXTENSION_LIMIT_SHORT):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
                consecutive_mismatches = 0
            else:
                mismatches += 1
                consecutive_mismatches += 1
                
                # Stop on too many consecutive mismatches or high error rate
                if (consecutive_mismatches >= 2 or 
                    (left_ext > config.T2_MISMATCH_THRESHOLD_EARLY and 
                     mismatches / left_ext > config.T2_MISMATCH_RATE)):
                    break
                left_ext += 1
        
        # Enhanced right extension
        right_ext = 0
        mismatches = 0
        consecutive_mismatches = 0
        
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < config.T2_EXTENSION_LIMIT_SHORT):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
                consecutive_mismatches = 0
            else:
                mismatches += 1
                consecutive_mismatches += 1
                
                if (consecutive_mismatches >= 2 or
                    (right_ext > config.T2_MISMATCH_THRESHOLD_EARLY and 
                     mismatches / right_ext > config.T2_MISMATCH_RATE)):
                    break
                right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    
    else:  # Reverse complement - exact matching only
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               left_ext < config.T2_EXTENSION_LIMIT_RC and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
        
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               right_ext < config.T2_EXTENSION_LIMIT_RC and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def extend_anchor_conservative_long(query, ref, qpos, rpos, strand, initial_len):
    """Conservative extension for long sequences with micro-optimizations"""
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1:
        # Slightly more aggressive extension for long sequences
        left_ext = 0
        mismatches = 0
        
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and 
               left_ext < config.T1_EXTENSION_LIMIT_LONG):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
            else:
                mismatches += 1
                # Slightly more tolerant stopping condition
                if (left_ext > config.T1_MISMATCH_THRESHOLD_EARLY and 
                    mismatches / left_ext > config.T1_MISMATCH_RATE * 1.1):  # 10% more tolerant
                    break
                left_ext += 1
        
        right_ext = 0
        mismatches = 0
        
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < config.T1_EXTENSION_LIMIT_LONG):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
            else:
                mismatches += 1
                if (right_ext > config.T1_MISMATCH_THRESHOLD_EARLY and 
                    mismatches / right_ext > config.T1_MISMATCH_RATE * 1.1):
                    break
                right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    
    else:  # Reverse complement
        left_ext = 0
        while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
               left_ext < config.T1_EXTENSION_LIMIT_RC and
               query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
            left_ext += 1
        
        right_ext = 0
        while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
               right_ext < config.T1_EXTENSION_LIMIT_RC and
               query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
            right_ext += 1
        
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def select_segments_conservative_short(segments, query, ref):
    """Conservative segment selection for short sequences with optimizations"""
    if not segments:
        return []
    
    # Enhanced validation and scoring
    valid_segments = []
    for seg in list(set(segments)):
        qstart, qend, rstart, rend = seg
        if (qend - qstart >= config.MIN_SEGMENT_LENGTH and 
            qend <= len(query) and rend <= len(ref) and
            qstart >= 0 and rstart >= 0 and qend > qstart):
            
            # Pre-validate with stricter criteria for short sequences
            if is_segment_valid(query, ref, qstart, qend, rstart, rend):
                # Additional quality score
                segment_score = qend - qstart
                valid_segments.append((seg, segment_score))
    
    if not valid_segments:
        return []
    
    # Sort by position, then by score
    valid_segments.sort(key=lambda x: (x[0][0], -x[1]))
    segments_only = [seg for seg, score in valid_segments]
    
    # Greedy selection
    selected = []
    last_query_end = 0
    
    for segment in segments_only:
        qstart, qend, rstart, rend = segment
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    
    return conservative_gap_fill_short(selected, query, ref)

def select_segments_long(segments, query, ref):
    """Segment selection for long sequences with optimizations"""
    valid_segments = list(set(segments))
    valid_segments.sort(key=lambda x: x[0])
    
    # Enhanced validation
    validated = []
    for seg in valid_segments:
        qstart, qend, rstart, rend = seg
        if (qend - qstart >= config.MIN_SEGMENT_LENGTH and 
            qend <= len(query) and rend <= len(ref) and
            qstart >= 0 and rstart >= 0 and qend > qstart and
            is_segment_valid(query, ref, qstart, qend, rstart, rend)):
            validated.append(seg)
    
    if not validated:
        return []
    
    # Optimized greedy selection
    selected = []
    last_query_end = 0
    
    for segment in validated:
        qstart, qend, rstart, rend = segment
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    
    return conservative_gap_fill_long(selected, query, ref)

def conservative_gap_fill_short(segments, query, ref):
    """Enhanced gap filling for short sequences"""
    if not segments:
        return []
    
    filled = []
    qlen = len(query)
    
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        
        # More aggressive start extension for short sequences
        if i == 0 and qstart > 0:
            gap_size = qstart
            if gap_size <= config.T2_GAP_FILL_START_MAX and rstart >= gap_size:
                # Allow up to 1 mismatch for small gaps
                mismatches = sum(1 for j in range(gap_size) 
                               if query[j] != ref[rstart - gap_size + j])
                
                if mismatches <= max(1, gap_size * 0.1):  # Allow 10% mismatches
                    extended_segment = (0, qend, rstart - gap_size, rend)
                    if is_segment_valid(query, ref, *extended_segment):
                        qstart, rstart = 0, rstart - gap_size
        
        # Enhanced gap filling between segments
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            
            if 0 < gap_size <= config.T2_GAP_FILL_BETWEEN_MAX:
                prev_segment = filled[-1]
                mismatches = 0
                valid_extension = True
                
                for j in range(gap_size):
                    if (prev_qend + j >= len(query) or 
                        prev_segment[3] + j >= len(ref)):
                        valid_extension = False
                        break
                    if query[prev_qend + j] != ref[prev_segment[3] + j]:
                        mismatches += 1
                
                # Allow small number of mismatches
                if (valid_extension and 
                    mismatches <= max(1, gap_size * 0.15)):
                    merged_segment = (prev_segment[0], qend, prev_segment[2], rend)
                    if is_segment_valid(query, ref, *merged_segment):
                        filled[-1] = merged_segment
                        continue
        
        filled.append((qstart, qend, rstart, rend))
    
    # Enhanced end extension
    if filled and filled[-1][1] < qlen:
        last_segment = filled[-1]
        gap_size = qlen - last_segment[1]
        if gap_size <= config.T2_GAP_FILL_END_MAX and last_segment[3] + gap_size <= len(ref):
            mismatches = sum(1 for j in range(gap_size)
                           if query[last_segment[1] + j] != ref[last_segment[3] + j])
            
            if mismatches <= max(1, gap_size * 0.1):
                extended_segment = (last_segment[0], qlen, last_segment[2], last_segment[3] + gap_size)
                if is_segment_valid(query, ref, *extended_segment):
                    filled[-1] = extended_segment
    
    return validate_no_overlaps(filled)

def conservative_gap_fill_long(segments, query, ref):
    """Enhanced gap filling for long sequences"""
    if not segments:
        return []
    
    filled = []
    qlen = len(query)
    
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        
        # Enhanced start extension
        if i == 0 and qstart > 0:
            gap_size = qstart
            if gap_size <= config.T1_GAP_FILL_START_MAX and rstart >= gap_size:
                mismatches = sum(1 for j in range(gap_size) 
                               if query[j] != ref[rstart - gap_size + j])
                
                # Slightly more tolerant for long sequences
                if mismatches <= max(1, gap_size * 0.08):
                    extended_segment = (0, qend, rstart - gap_size, rend)
                    if is_segment_valid(query, ref, *extended_segment):
                        qstart, rstart = 0, rstart - gap_size
        
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            
            if 0 < gap_size <= config.T1_GAP_FILL_BETWEEN_MAX:
                prev_segment = filled[-1]
                mismatches = 0
                valid_extension = True
                
                for j in range(gap_size):
                    if (prev_qend + j >= len(query) or 
                        prev_segment[3] + j >= len(ref)):
                        valid_extension = False
                        break
                    if query[prev_qend + j] != ref[prev_segment[3] + j]:
                        mismatches += 1
                
                if (valid_extension and 
                    mismatches <= max(1, gap_size * 0.1)):
                    merged_segment = (prev_segment[0], qend, prev_segment[2], rend)
                    if is_segment_valid(query, ref, *merged_segment):
                        filled[-1] = merged_segment
                        continue
        
        filled.append((qstart, qend, rstart, rend))
    
    # Enhanced end extension
    if filled and filled[-1][1] < qlen:
        last_segment = filled[-1]
        gap_size = qlen - last_segment[1]
        if gap_size <= config.T1_GAP_FILL_END_MAX and last_segment[3] + gap_size <= len(ref):
            mismatches = sum(1 for j in range(gap_size)
                           if query[last_segment[1] + j] != ref[last_segment[3] + j])
            
            if mismatches <= max(1, gap_size * 0.08):
                extended_segment = (last_segment[0], qlen, last_segment[2], last_segment[3] + gap_size)
                if is_segment_valid(query, ref, *extended_segment):
                    filled[-1] = extended_segment
    
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
