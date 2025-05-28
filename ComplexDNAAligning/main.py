from evaluate import calculate_value, seq2hashtable_multi_test, calculate_distance
import config
import time

def find_alignment_unified(query, ref, kmer_sizes, params):
    all_anchors = []
    if len(query) < config.SHORT_SEQUENCE_THRESHOLD:
        shift = 1 
        for kmer_size in kmer_sizes:
            if kmer_size <= 5:
                anchors = seq2hashtable_multi_test(ref, query, kmer_size, config.T2_SMALL_KMER_SHIFT)
                all_anchors.extend(anchors)
            elif kmer_size <= 8:
                anchors = seq2hashtable_multi_test(ref, query, kmer_size, config.T2_MEDIUM_KMER_SHIFT)
                all_anchors.extend(anchors)
            else:
                anchors = seq2hashtable_multi_test(ref, query, kmer_size, shift)
                all_anchors.extend(anchors)
    else:
        shift = config.DENSE_SHIFT
        for kmer_size in kmer_sizes:
            anchors = seq2hashtable_multi_test(ref, query, kmer_size, shift)
            all_anchors.extend(anchors)
    if not all_anchors:
        return []
    
    segments = []
    for query_pos, ref_pos, strand, match_len in all_anchors:
        segment = extend_anchor_unified(query, ref, query_pos, ref_pos, strand, match_len, params)
        if segment and segment[1] - segment[0] >= config.MIN_SEGMENT_LENGTH:
            if is_segment_valid(query, ref, *segment):
                segments.append(segment)
    
    return select_segments_unified(segments, query, ref, params)

def extend_anchor_unified(query, ref, qpos, rpos, strand, initial_len, params):
    qlen, rlen = len(query), len(ref)
    rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    if strand == 1: 
        left_ext = 0
        mismatches = 0
        while (qpos - left_ext > 0 and rpos - left_ext > 0 and left_ext < params['extension_limit']):
            if query[qpos - left_ext - 1] == ref[rpos - left_ext - 1]:
                left_ext += 1
            else:
                mismatches += 1
                if left_ext > params['mismatch_threshold'] and mismatches / left_ext > params['mismatch_rate']:
                    break
                left_ext += 1
        right_ext = 0
        mismatches = 0
        while (qpos + initial_len + right_ext < qlen and 
               rpos + initial_len + right_ext < rlen and 
               right_ext < params['extension_limit']):
            if query[qpos + initial_len + right_ext] == ref[rpos + initial_len + right_ext]:
                right_ext += 1
            else:
                mismatches += 1
                if right_ext > params['mismatch_threshold'] and mismatches / right_ext > params['mismatch_rate']:
                    break
                right_ext += 1
        if qlen < config.SHORT_SEQUENCE_THRESHOLD:
            if qpos - left_ext <= config.T1_BOUNDARY_TOLERANCE and rpos - left_ext >= 0:
                boundary_left = min(qpos, rpos)
                left_ext = max(left_ext, boundary_left)
            if qpos + initial_len + right_ext >= qlen - config.T1_BOUNDARY_TOLERANCE and rpos + initial_len + right_ext + (qlen - qpos - initial_len - right_ext) <= rlen:
                boundary_right = qlen - qpos - initial_len
                right_ext = max(right_ext, boundary_right)
        return (qpos - left_ext, qpos + initial_len + right_ext, 
                rpos - left_ext, rpos + initial_len + right_ext)
    else: 
        if qlen < config.SHORT_SEQUENCE_THRESHOLD:
            left_ext = 0
            left_mismatches = 0
            while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
                   left_ext < params['extension_limit_rc']):
                expected_char = rc_map.get(ref[rpos + initial_len + left_ext], 'N')
                if query[qpos - left_ext - 1] == expected_char:
                    left_ext += 1
                else:
                    left_mismatches += 1
                    if left_ext > config.T1_RC_MIN_LENGTH and left_mismatches / left_ext > config.T1_RC_MISMATCH_RATE:
                        break
                    left_ext += 1
            right_ext = 0
            right_mismatches = 0
            while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
                   right_ext < params['extension_limit_rc']):
                expected_char = rc_map.get(ref[rpos - right_ext - 1], 'N')
                if query[qpos + initial_len + right_ext] == expected_char:
                    right_ext += 1
                else:
                    right_mismatches += 1
                    if right_ext > config.T1_RC_MIN_LENGTH and right_mismatches / right_ext > config.T1_RC_MISMATCH_RATE:
                        break
                    right_ext += 1
        else:
            left_ext = 0
            while (qpos - left_ext > 0 and rpos + initial_len + left_ext < rlen and
                   left_ext < params['extension_limit_rc'] and
                   query[qpos - left_ext - 1] == rc_map.get(ref[rpos + initial_len + left_ext], 'N')):
                left_ext += 1
            
            right_ext = 0
            while (qpos + initial_len + right_ext < qlen and rpos - right_ext > 0 and
                   right_ext < params['extension_limit_rc'] and
                   query[qpos + initial_len + right_ext] == rc_map.get(ref[rpos - right_ext - 1], 'N')):
                right_ext += 1
        return (qpos - left_ext, qpos + initial_len + right_ext,
                rpos - right_ext, rpos + initial_len + left_ext)

def select_segments_unified(segments, query, ref, params):
    if not segments:
        return []
    valid_segments = list(set(segments))
    valid_segments.sort(key=lambda x: x[0])
    selected = []
    last_query_end = 0
    for segment in valid_segments:
        qstart, qend, rstart, rend = segment
        if qstart >= last_query_end:
            selected.append(segment)
            last_query_end = qend
    return gap_fill_unified(selected, query, ref, params)

def gap_fill_unified(segments, query, ref, params):
    if not segments:
        return []
    filled = []
    qlen = len(query)
    for i, segment in enumerate(segments):
        qstart, qend, rstart, rend = segment
        if i == 0 and qstart > 0:
            gap_size = qstart
            if gap_size <= params['gap_fill_start_max'] and rstart >= gap_size:
                can_extend = all(query[j] == ref[rstart - gap_size + j] 
                               for j in range(gap_size))
                if can_extend:
                    extended_segment = (0, qend, rstart - gap_size, rend)
                    if is_segment_valid(query, ref, *extended_segment):
                        qstart, rstart = 0, rstart - gap_size
        if filled:
            prev_qend = filled[-1][1]
            gap_size = qstart - prev_qend
            if 0 < gap_size <= params['gap_fill_between_max']:
                prev_segment = filled[-1]
                can_extend = all(
                    prev_qend + j < len(query) and 
                    prev_segment[3] + j < len(ref) and
                    query[prev_qend + j] == ref[prev_segment[3] + j]
                    for j in range(gap_size)
                )
                if can_extend:
                    merged_segment = (prev_segment[0], qend, prev_segment[2], rend)
                    if is_segment_valid(query, ref, *merged_segment):
                        filled[-1] = merged_segment
                        continue
        filled.append((qstart, qend, rstart, rend))
   
    if filled and filled[-1][1] < qlen:
        last_segment = filled[-1]
        gap_size = qlen - last_segment[1]
        if gap_size <= params['gap_fill_end_max'] and last_segment[3] + gap_size <= len(ref):
            can_extend = all(query[last_segment[1] + j] == ref[last_segment[3] + j]
                           for j in range(gap_size))
            if can_extend:
                extended_segment = (last_segment[0], qlen, last_segment[2], last_segment[3] + gap_size)
                if is_segment_valid(query, ref, *extended_segment):
                    filled[-1] = extended_segment
   
    if len(query) < config.SHORT_SEQUENCE_THRESHOLD:
        final_filled = []
        for i, segment in enumerate(filled):
            qstart, qend, rstart, rend = segment
            if final_filled:
                prev_qend = final_filled[-1][1]
                gap_size = qstart - prev_qend
                if 0 < gap_size <= config.T2_SMALL_GAP_MAX:
                    prev_segment = final_filled[-1]
                    mismatches = sum(1 for j in range(gap_size) 
                                   if prev_qend + j < len(query) and 
                                      prev_segment[3] + j < len(ref) and
                                      query[prev_qend + j] != ref[prev_segment[3] + j])
                    
                    if mismatches / gap_size <= config.T2_SMALL_GAP_MISMATCH_RATE: 
                        merged_segment = (prev_segment[0], qend, prev_segment[2], rend)
                        if is_segment_valid(query, ref, *merged_segment):
                            final_filled[-1] = merged_segment
                            continue
            final_filled.append(segment)
        return validate_no_overlaps(final_filled)
    else:
        final_filled = []
        for i, segment in enumerate(filled):
            qstart, qend, rstart, rend = segment
            
            if final_filled:
                prev_qend = final_filled[-1][1]
                gap_size = qstart - prev_qend
                if 0 < gap_size <= config.T1_LARGE_GAP_MAX: 
                    prev_segment = final_filled[-1]
                   
                    matches = sum(1 for j in range(gap_size) 
                                if prev_qend + j < len(query) and 
                                   prev_segment[3] + j < len(ref) and
                                   query[prev_qend + j] == ref[prev_segment[3] + j])
                    
                    if matches / gap_size >= config.T1_LARGE_GAP_MATCH_RATE: 
                        merged_segment = (prev_segment[0], qend, prev_segment[2], rend)
                       
                        if qend - prev_segment[0] >= config.MIN_SEGMENT_LENGTH:
                            try:
                                edit_dist = calculate_distance(ref, query, prev_segment[2], rend, prev_segment[0], qend)
                                if edit_dist / (qend - prev_segment[0]) <= config.T1_MERGED_SEGMENT_EDIT_RATE: 
                                    final_filled[-1] = merged_segment
                                    continue
                            except:
                                pass
            final_filled.append(segment)
        return validate_no_overlaps(final_filled)

def find_alignment_tuples(query, ref):
    if len(query) < config.SHORT_SEQUENCE_THRESHOLD:
        kmer_sizes = config.T2_KMER_SIZES
        params = {
            'extension_limit': config.T2_EXTENSION_LIMIT_SHORT,
            'extension_limit_rc': config.T2_EXTENSION_LIMIT_RC,
            'mismatch_threshold': config.T2_MISMATCH_THRESHOLD_EARLY,
            'mismatch_rate': config.T2_MISMATCH_RATE,
            'gap_fill_start_max': config.T2_GAP_FILL_START_MAX,
            'gap_fill_between_max': config.T2_GAP_FILL_BETWEEN_MAX,
            'gap_fill_end_max': config.T2_GAP_FILL_END_MAX
        }
    else:
        kmer_sizes = config.T1_KMER_SIZES
        params = {
            'extension_limit': config.T1_EXTENSION_LIMIT_LONG,
            'extension_limit_rc': config.T1_EXTENSION_LIMIT_RC,
            'mismatch_threshold': config.T1_MISMATCH_THRESHOLD_EARLY,
            'mismatch_rate': config.T1_MISMATCH_RATE,
            'gap_fill_start_max': config.T1_GAP_FILL_START_MAX,
            'gap_fill_between_max': config.T1_GAP_FILL_BETWEEN_MAX,
            'gap_fill_end_max': config.T1_GAP_FILL_END_MAX
        }
    
    return find_alignment_unified(query, ref, kmer_sizes, params)

def is_segment_valid(query, ref, qstart, qend, rstart, rend):
    if qend - qstart < config.MIN_SEGMENT_LENGTH:
        return False
    try:
        edit_dist = calculate_distance(ref, query, rstart, rend, qstart, qend)
        edit_rate = edit_dist / (qend - qstart)
        return edit_rate <= config.MAX_EDIT_RATE
    except:
        return False

def validate_no_overlaps(segments):
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
    start_time = time.time()
    
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
    
    print(f"Data loading completed in {time.time() - start_time:.3f} seconds")
    
    t1_start = time.time()
    t1 = find_alignment_tuples(q1, r1)
    t1_time = time.time() - t1_start
   
    t2_start = time.time()
    t2 = find_alignment_tuples(q2, r2)
    t2_time = time.time() - t2_start
    
    t1_score = calculate_value(str(t1), r1, q1)
    t2_score = calculate_value(str(t2), r2, q2)
    
    print('t1: ', t1)
    print('t2: ', t2)
    print(f't1 score: {t1_score}')
    print(f't2 score: {t2_score}')
    print(f'T1 alignment time: {t1_time:.3f} seconds')
    print(f'T2 alignment time: {t2_time:.3f} seconds')
    
if __name__ == "__main__":
    main()
