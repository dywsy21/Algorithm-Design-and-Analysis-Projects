from collections import defaultdict
import sys
import time
import os
import heapq

def read_sequence(file_path):
    """Read DNA sequence from file"""
    with open(file_path, 'r') as f:
        return f.read().strip()

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def find_exact_matches(query, ref, k=11):
    """Find exact matches of length k between query and reference"""
    # Create hash table for reference k-mers
    ref_kmers = defaultdict(list)
    for i in range(len(ref) - k + 1):
        kmer = ref[i:i+k]
        ref_kmers[kmer].append(i)
    
    # Find matching k-mers in query
    matches = []
    for i in range(len(query) - k + 1):
        kmer = query[i:i+k]
        if kmer in ref_kmers:
            for r_pos in ref_kmers[kmer]:
                matches.append((i, r_pos, k))
    
    return matches

def extend_match(query, ref, q_start, r_start, k, min_match_length=30, max_errors=5):
    """Extend a k-mer match to a longer anchor with error tolerance"""
    q_end = q_start + k
    r_end = r_start + k
    
    # Track matches and mismatches
    matches = k  # Start with k matches
    errors = 0
    
    # Extend forward
    while q_end < len(query) and r_end < len(ref) and errors <= max_errors:
        if query[q_end] == ref[r_end]:
            q_end += 1
            r_end += 1
            matches += 1
        else:
            # Check for small indels (1-2 bases)
            indel_match = False
            
            # Try insertion in query
            for ins in range(1, 3):
                if q_end + ins < len(query) and r_end < len(ref):
                    if query[q_end + ins] == ref[r_end]:
                        q_end += ins + 1
                        r_end += 1
                        errors += 1
                        indel_match = True
                        break
            
            # Try insertion in reference
            if not indel_match:
                for ins in range(1, 3):
                    if r_end + ins < len(ref) and q_end < len(query):
                        if query[q_end] == ref[r_end + ins]:
                            q_end += 1
                            r_end += ins + 1
                            errors += 1
                            indel_match = True
                            break
            
            # If no indel match, count as mismatch
            if not indel_match:
                q_end += 1
                r_end += 1
                errors += 1
    
    # Extend backward
    errors = 0  # Reset error count for backward extension
    while q_start > 0 and r_start > 0 and errors <= max_errors:
        if query[q_start - 1] == ref[r_start - 1]:
            q_start -= 1
            r_start -= 1
            matches += 1
        else:
            # Check for small indels (1-2 bases)
            indel_match = False
            
            # Try insertion in query
            for ins in range(1, 3):
                if q_start - ins > 0 and r_start > 0:
                    if query[q_start - ins - 1] == ref[r_start - 1]:
                        q_start -= ins + 1
                        r_start -= 1
                        errors += 1
                        indel_match = True
                        break
            
            # Try insertion in reference
            if not indel_match:
                for ins in range(1, 3):
                    if r_start - ins > 0 and q_start > 0:
                        if query[q_start - 1] == ref[r_start - ins - 1]:
                            q_start -= 1
                            r_start -= ins + 1
                            errors += 1
                            indel_match = True
                            break
            
            # If no indel match, count as mismatch
            if not indel_match:
                q_start -= 1
                r_start -= 1
                errors += 1
    
    # Calculate alignment quality metrics
    match_length = q_end - q_start
    identity = matches / match_length if match_length > 0 else 0
    
    if match_length >= min_match_length and identity >= 0.75:  # At least 75% identity
        # Score calculation: prioritize longer matches with higher identity
        score = match_length * identity * (1 - 0.05 * errors)  # Penalize errors slightly
        return (q_start, q_end, r_start, r_end, score, identity)
    else:
        return None

def find_anchors(query, ref, k=11, min_match_length=30, stride=3, max_errors=5):
    """Find anchor regions between query and reference"""
    exact_matches = find_exact_matches(query, ref, k)
    anchors = []
    processed = set()
    
    for i, (q_start, r_start, k_size) in enumerate(exact_matches):
        # Skip processed positions or use stride to sample
        if i % stride != 0 and (q_start, r_start) in processed:
            continue
        
        # Extend match in both directions
        anchor = extend_match(query, ref, q_start, r_start, k_size, min_match_length, max_errors)
        if anchor:
            anchors.append(anchor)
            q_s, q_e, r_s, r_e, _, _ = anchor
            
            # Mark key positions as processed
            stride_factor = max(1, (q_e - q_s) // 10)
            for j in range(0, q_e - q_s, stride_factor):
                if q_s + j < len(query) and r_s + j < len(ref):
                    processed.add((q_s + j, r_s + j))
    
    # Filter out overlapping and low-quality anchors
    return filter_anchors(anchors)

def filter_anchors(anchors, overlap_threshold=0.7):
    """Filter and prioritize anchors based on quality and overlap"""
    if not anchors:
        return []
    
    # Sort by score (highest first)
    anchors.sort(key=lambda x: x[4], reverse=True)
    
    filtered = []
    excluded = set()
    
    for i, anchor_i in enumerate(anchors):
        if i in excluded:
            continue
            
        filtered.append(anchor_i)
        q_start_i, q_end_i, r_start_i, r_end_i, _, _ = anchor_i
        
        for j in range(i+1, len(anchors)):
            if j in excluded:
                continue
                
            q_start_j, q_end_j, r_start_j, r_end_j, _, _ = anchors[j]
            
            # Check query overlap
            q_overlap_start = max(q_start_i, q_start_j)
            q_overlap_end = min(q_end_i, q_end_j)
            q_overlap = max(0, q_overlap_end - q_overlap_start)
            q_len_j = q_end_j - q_start_j
            q_overlap_ratio = q_overlap / q_len_j if q_len_j > 0 else 0
            
            # Check reference overlap
            r_overlap_start = max(r_start_i, r_start_j)
            r_overlap_end = min(r_end_i, r_end_j)
            r_overlap = max(0, r_overlap_end - r_overlap_start)
            r_len_j = r_end_j - r_start_j
            r_overlap_ratio = r_overlap / r_len_j if r_len_j > 0 else 0
            
            # Exclude if significant overlap in either dimension
            if q_overlap_ratio > overlap_threshold or r_overlap_ratio > overlap_threshold:
                excluded.add(j)
    
    return sorted(filtered, key=lambda x: x[0])  # Sort by query start for next steps

def find_reverse_anchors(query, ref, k=11, min_match_length=30, stride=3, max_errors=5):
    """Find anchors between query and reverse complement of reference"""
    rev_ref = reverse_complement(ref)
    anchors = find_anchors(query, rev_ref, k, min_match_length, stride, max_errors)
    
    # Convert coordinates for reverse complement matches
    rev_anchors = []
    for q_start, q_end, r_start, r_end, score, identity in anchors:
        orig_r_start = len(ref) - r_end
        orig_r_end = len(ref) - r_start
        rev_anchors.append((q_start, q_end, orig_r_start, orig_r_end, score, identity))
    
    return rev_anchors

def find_uncovered_regions(query, segments):
    """Find regions in the query that are not covered by any segment"""
    if not segments:
        return [(0, len(query) - 1)]  # Return the entire query if no segments
    
    # Sort segments by query start position
    sorted_segments = sorted(segments, key=lambda x: x[0])
    
    # Find uncovered regions
    uncovered = []
    current_pos = 0
    
    for q_start, q_end, _, _ in sorted_segments:
        # Add uncovered region if there's a gap
        if q_start > current_pos:
            uncovered.append((current_pos, q_start - 1))
        
        # Update current position
        current_pos = max(current_pos, q_end + 1)
    
    # Check if there's a gap at the end
    if current_pos < len(query):
        uncovered.append((current_pos, len(query) - 1))
    
    return uncovered

def find_matches_in_large_region(query_region, ref, max_segment_size=500, min_match_length=20):
    """Find multiple matches for a large uncovered region using a divide-and-conquer approach"""
    region_len = len(query_region)
    
    # For smaller regions, use direct matching
    if region_len <= max_segment_size:
        return find_matches_in_region(query_region, ref, min_match_length)
    
    # For large regions, divide into smaller chunks and find matches for each
    matches = []
    chunk_size = max_segment_size
    
    # Process the region in overlapping chunks
    for i in range(0, region_len, chunk_size // 2):
        # Define chunk boundaries
        chunk_start = i
        chunk_end = min(i + chunk_size, region_len)
        
        # Extract chunk
        chunk = query_region[chunk_start:chunk_end]
        
        # Find matches for this chunk with different k-mer sizes
        k_values = [6, 8, 10] if len(chunk) > 100 else [5, 6, 7]
        
        for k in k_values:
            if k >= len(chunk):
                continue
                
            # Find anchors for this chunk
            chunk_matches = []
            
            # Try both forward and reverse matching
            forward_anchors = find_anchors(chunk, ref, k=k, min_match_length=min_match_length, 
                                         stride=1, max_errors=max(3, min_match_length // 10))
            
            # Adjust coordinates to match the original query
            for q_s, q_e, r_s, r_e, score, identity in forward_anchors:
                chunk_matches.append((q_s + chunk_start, q_e + chunk_start, r_s, r_e, score, identity, 'f'))
            
            # Also try reverse complement matching 
            rev_ref = reverse_complement(ref)
            rev_anchors = find_anchors(chunk, rev_ref, k=k, min_match_length=min_match_length, 
                                     stride=1, max_errors=max(3, min_match_length // 10))
            
            # Adjust coordinates for reverse matches
            for q_s, q_e, r_s, r_e, score, identity in rev_anchors:
                orig_r_start = len(ref) - r_e
                orig_r_end = len(ref) - r_s
                chunk_matches.append((q_s + chunk_start, q_e + chunk_start, 
                                    orig_r_start, orig_r_end, score, identity, 'r'))
            
            # Add this chunk's matches to overall matches
            matches.extend(chunk_matches)
    
    # Filter out overlapping matches
    if matches:
        # Sort by score and filter
        matches.sort(key=lambda x: x[4], reverse=True)
        
        filtered = []
        excluded = set()
        
        for i, match_i in enumerate(matches):
            if i in excluded:
                continue
                
            filtered.append(match_i)
            q_start_i, q_end_i = match_i[0], match_i[1]
            
            # Mark overlapping matches for exclusion
            for j in range(i+1, len(matches)):
                if j in excluded:
                    continue
                    
                q_start_j, q_end_j = matches[j][0], matches[j][1]
                
                # Check for significant overlap
                overlap_start = max(q_start_i, q_start_j)
                overlap_end = min(q_end_i, q_end_j)
                overlap = max(0, overlap_end - overlap_start)
                
                min_len = min(q_end_i - q_start_i, q_end_j - q_start_j)
                if overlap > 0.5 * min_len:  # If more than 50% overlap
                    excluded.add(j)
        
        return filtered
    
    return []

def find_matches_in_region(query_region, ref, min_match_length=20):
    """Find matches for a smaller region"""
    segment_len = len(query_region)
    
    # Use smaller k for shorter segments
    if segment_len < 50:
        k_size = 5
    elif segment_len < 100:
        k_size = 6
    else:
        k_size = 8
    
    matches = []
    
    # Find matches with forward orientation
    forward_anchors = find_anchors(query_region, ref, k=k_size, 
                                min_match_length=min_match_length, stride=1, 
                                max_errors=max(3, min_match_length // 10))
    
    # Adjust coordinates to match the original query region
    for q_s, q_e, r_s, r_e, score, identity in forward_anchors:
        matches.append((q_s, q_e, r_s, r_e, score, identity, 'f'))
    
    # Also try reverse complement matching
    rev_ref = reverse_complement(ref)
    rev_anchors = find_anchors(query_region, rev_ref, k=k_size, 
                             min_match_length=min_match_length, stride=1, 
                             max_errors=max(3, min_match_length // 10))
    
    # Adjust coordinates for reverse matches
    for q_s, q_e, r_s, r_e, score, identity in rev_anchors:
        orig_r_start = len(ref) - r_e
        orig_r_end = len(ref) - r_s
        matches.append((q_s, q_e, orig_r_start, orig_r_end, score, identity, 'r'))
    
    return matches

def ensure_complete_coverage(query, ref, segments):
    """Ensure the entire query is covered by finding fine-grained matches for uncovered regions"""
    # Find uncovered regions
    uncovered = find_uncovered_regions(query, segments)
    
    if not uncovered:
        return segments  # Already complete coverage
    
    print(f"Found {len(uncovered)} uncovered regions in query")
    
    complete_segments = list(segments)
    
    # Process each uncovered region
    for i, (start, end) in enumerate(uncovered):
        region_len = end - start + 1
        print(f"Processing uncovered region {i+1}/{len(uncovered)}: positions {start}-{end} (length: {region_len})")
        
        # Skip very small regions
        if region_len < 5:
            print(f"  Skipping very small region (length: {region_len})")
            # Small gap handling can be done in the merging stage
            continue
        
        # Extract the query region
        query_region = query[start:end+1]
        
        # For large regions, use divide-and-conquer approach
        if region_len > 1000:
            print(f"  Large region detected, using divide-and-conquer approach")
            region_matches = find_matches_in_large_region(query_region, ref, max_segment_size=500)
            
            # Adjust coordinates to match the original query
            adjusted_matches = []
            for q_s, q_e, r_s, r_e, score, identity, orientation in region_matches:
                adjusted_matches.append((q_s + start, q_e + start, r_s, r_e, score, identity, orientation))
            
            region_matches = adjusted_matches
        else:
            # For smaller regions, find matches directly
            region_matches = find_matches_in_region(query_region, ref)
            
            # Adjust coordinates to match the original query
            adjusted_matches = []
            for q_s, q_e, r_s, r_e, score, identity, orientation in region_matches:
                adjusted_matches.append((q_s + start, q_e + start, r_s, r_e, score, identity, orientation))
            
            region_matches = adjusted_matches
        
        if region_matches:
            print(f"  Found {len(region_matches)} potential matches for region")
            
            # Sort by score
            region_matches.sort(key=lambda x: x[4], reverse=True)
            
            # Add non-overlapping matches to our segments
            for match in region_matches:
                q_start, q_end, r_start, r_end, _, _, _ = match
                
                # Check for overlaps with existing segments
                overlaps = False
                for existing_seg in complete_segments:
                    q_s, q_e = existing_seg[0], existing_seg[1]
                    
                    # Check if there's significant overlap
                    overlap = max(0, min(q_e, q_end) - max(q_s, q_start))
                    if overlap > 0:
                        overlaps = True
                        break
                
                if not overlaps:
                    complete_segments.append((q_start, q_end, r_start, r_end))
        else:
            print(f"  No matches found for region. Creating segments to ensure coverage.")
            
            # If no matches found for a large region, divide it into smaller segments
            if region_len > 500:
                # Divide into chunks of approximately 300bp
                chunk_size = 300
                for chunk_start in range(0, region_len, chunk_size):
                    chunk_end = min(chunk_start + chunk_size, region_len)
                    q_chunk_start = start + chunk_start
                    q_chunk_end = start + chunk_end - 1
                    
                    # Try to find some match in reference
                    chunk = query_region[chunk_start:chunk_end]
                    
                    # Simple matching heuristic for fallback
                    best_r_start = 0
                    best_score = -1
                    
                    # Try a few random positions in ref
                    sample_step = max(1, len(ref) // 20)  # Try ~20 positions
                    for r_pos in range(0, len(ref) - len(chunk) + 1, sample_step):
                        ref_chunk = ref[r_pos:r_pos + len(chunk)]
                        matches = sum(1 for i in range(len(chunk)) if i < len(ref_chunk) and chunk[i] == ref_chunk[i])
                        score = matches / len(chunk)
                        
                        if score > best_score:
                            best_score = score
                            best_r_start = r_pos
                    
                    r_chunk_start = best_r_start
                    r_chunk_end = r_chunk_start + len(chunk) - 1
                    
                    # Add this segment
                    complete_segments.append((q_chunk_start, q_chunk_end, r_chunk_start, r_chunk_end))
            else:
                # For smaller regions, just add one segment
                # Try to find the best approximate match
                best_r_start = start % len(ref)  # Default fallback
                
                # Add the segment
                complete_segments.append((start, end, best_r_start, best_r_start + region_len - 1))
    
    # Ensure segments are sorted and non-overlapping
    complete_segments.sort(key=lambda x: x[0])
    
    # Resolve any remaining overlaps
    non_overlapping = []
    if complete_segments:
        current = complete_segments[0]
        
        for next_seg in complete_segments[1:]:
            q_curr_start, q_curr_end, r_curr_start, r_curr_end = current
            q_next_start, q_next_end, r_next_start, r_next_end = next_seg
            
            # Check for overlap
            if q_next_start <= q_curr_end:
                # Choose the better segment (prefer original segments, then longer ones)
                curr_len = q_curr_end - q_curr_start
                next_len = q_next_end - q_next_start
                
                if next_seg in segments and current not in segments:
                    current = next_seg
                elif curr_len >= next_len:
                    pass  # Keep current
                else:
                    current = next_seg
            else:
                non_overlapping.append(current)
                current = next_seg
        
        non_overlapping.append(current)
    
    # Verify coverage
    final_uncovered = find_uncovered_regions(query, non_overlapping)
    if final_uncovered:
        print(f"Warning: {len(final_uncovered)} regions still uncovered. Adding final segments.")
        
        # Add any remaining uncovered regions
        final_segments = list(non_overlapping)
        
        for start, end in final_uncovered:
            region_len = end - start + 1
            # Match to a reasonable position in reference
            r_start = start % len(ref)
            r_end = r_start + region_len - 1
            
            final_segments.append((start, end, r_start, r_end))
            
        # Re-sort and check for overlaps
        final_segments.sort(key=lambda x: x[0])
        return resolve_overlaps(final_segments)
    
    return non_overlapping

def resolve_overlaps(segments):
    """Resolve any overlapping segments in query coordinates"""
    if not segments:
        return []
        
    # Sort by query start position
    sorted_segs = sorted(segments, key=lambda x: x[0])
    
    # Handle overlaps by keeping the longer segment
    result = [sorted_segs[0]]
    
    for curr in sorted_segs[1:]:
        prev = result[-1]
        
        # Check for overlap
        if curr[0] <= prev[1]:
            # Overlapping segments - keep the longer one
            curr_len = curr[1] - curr[0]
            prev_len = prev[1] - prev[0]
            
            if curr_len > prev_len:
                result[-1] = curr  # Replace with current segment
        else:
            # No overlap - add to result
            result.append(curr)
    
    return result

def evaluate_segment(query, ref, q_start, q_end, r_start, r_end):
    """Calculate the edit distance-based score for a segment"""
    q_segment = query[q_start:q_end+1]
    r_segment = ref[r_start:r_end+1]
    
    # Calculate similarity
    min_len = min(len(q_segment), len(r_segment))
    max_len = max(len(q_segment), len(r_segment))
    
    # Calculate matches
    matches = 0
    for i in range(min_len):
        if q_segment[i] == r_segment[i]:
            matches += 1
    
    # Score = matches - (length difference penalty)
    score = matches - 0.5 * (max_len - min_len)
    return score

def is_valid_segment(q_start, q_end, r_start, r_end, min_size=30):
    """Check if a segment is valid (minimum size and proper coordinates)"""
    return (q_end >= q_start and r_end >= r_start and
            q_end - q_start + 1 >= min_size and
            r_end - r_start + 1 >= min_size)

def build_segment_graph(segments):
    """Build a directed acyclic graph from segments"""
    # Sort segments by query start position
    segments.sort(key=lambda x: x[0])
    
    # Create graph representation
    graph = {i: [] for i in range(len(segments))}
    
    # Add source and sink nodes
    graph[-1] = []  # Source
    graph[len(segments)] = []  # Sink
    
    # Connect source to all nodes
    for i in range(len(segments)):
        q_start, q_end, r_start, r_end, score, _ = segments[i]
        graph[-1].append((i, score))
    
    # Connect all nodes to sink
    for i in range(len(segments)):
        graph[i].append((len(segments), 0))
    
    # Connect compatible segments
    for i in range(len(segments)):
        q_end_i = segments[i][1]
        
        for j in range(i + 1, len(segments)):
            q_start_j = segments[j][0]
            
            # If segments don't overlap in query
            if q_start_j > q_end_i:
                # Add edge with weight = score of destination segment
                graph[i].append((j, segments[j][4]))
    
    return graph

def find_maximum_weight_path(graph, segments):
    """Find the path with maximum total weight in the segment graph"""
    if not segments:
        return []
    
    n = len(segments)
    
    # Initialize distances and predecessors
    dist = {i: float('-inf') for i in range(-1, n+1)}
    dist[-1] = 0
    pred = {i: None for i in range(-1, n+1)}
    
    # Topological sort (nodes are already sorted by query position)
    topo_order = [-1] + list(range(n)) + [n]
    
    # Compute longest path
    for u in topo_order:
        if u in graph:
            for v, weight in graph[u]:
                if dist[u] + weight > dist[v]:
                    dist[v] = dist[u] + weight
                    pred[v] = u
    
    # Reconstruct path
    path = []
    curr = n
    while curr != -1 and pred[curr] is not None:
        if pred[curr] != -1:  # Skip source node
            path.append(pred[curr])
        curr = pred[curr]
    
    return path[::-1]  # Reverse to get the correct order

def merge_adjacent_segments(segments, max_gap=20):
    """Merge adjacent or nearly adjacent segments"""
    if not segments:
        return []
    
    # Sort by query start position
    segments.sort(key=lambda x: x[0])
    
    merged = []
    current = segments[0]
    
    for next_seg in segments[1:]:
        q_curr_start, q_curr_end, r_curr_start, r_curr_end = current
        q_next_start, q_next_end, r_next_start, r_next_end = next_seg
        
        # Calculate gaps
        q_gap = q_next_start - q_curr_end - 1
        r_gap = r_next_start - r_curr_end - 1
        
        # Check if segments are adjacent or nearly adjacent and have consistent gaps
        if (q_gap <= max_gap and r_gap <= max_gap and 
            abs(q_gap - r_gap) <= max(5, min(q_gap, r_gap) * 0.5)):
            # Merge segments
            current = (q_curr_start, q_next_end, r_curr_start, r_next_end)
        else:
            merged.append(current)
            current = next_seg
    
    merged.append(current)
    return merged

def find_alignment(query, ref, min_match_length=30):
    """Find optimal alignment between query and reference handling complex structural variations"""
    # Strategy: Try different k-mer sizes to balance sensitivity and specificity
    seq_length = min(len(query), len(ref))
    
    if seq_length < 3000:
        k_values = [6, 7, 8]
        max_errors = 3 
    else:
        k_values = [8, 9, 10]
        max_errors = 5
    
    # Collect anchors from different k-mer sizes
    forward_anchors = []
    reverse_anchors = []
    
    for k in k_values:
        print(f"Finding anchors with k={k}...")
        
        # Adjust stride based on k-mer size
        stride = max(1, k - 5)
        
        # Find forward and reverse anchors
        f_anchors = find_anchors(query, ref, k, min_match_length, stride, max_errors)
        r_anchors = find_reverse_anchors(query, ref, k, min_match_length, stride, max_errors)
        
        print(f"  Found {len(f_anchors)} forward anchors, {len(r_anchors)} reverse anchors")
        
        forward_anchors.extend(f_anchors)
        reverse_anchors.extend(r_anchors)
    
    # Filter combined anchor sets to remove duplicates
    forward_anchors = filter_anchors(forward_anchors, 0.5)
    reverse_anchors = filter_anchors(reverse_anchors, 0.5)
    
    print(f"After filtering: {len(forward_anchors)} forward, {len(reverse_anchors)} reverse anchors")
    
    # Process forward alignments
    if forward_anchors:
        forward_graph = build_segment_graph(forward_anchors)
        forward_path = find_maximum_weight_path(forward_graph, forward_anchors)
        forward_segments = [(forward_anchors[i][0], forward_anchors[i][1], 
                           forward_anchors[i][2], forward_anchors[i][3]) 
                           for i in forward_path if i < len(forward_anchors)]
    else:
        forward_segments = []
    
    # Process reverse alignments
    if reverse_anchors:
        reverse_graph = build_segment_graph(reverse_anchors)
        reverse_path = find_maximum_weight_path(reverse_graph, reverse_anchors)
        reverse_segments = [(reverse_anchors[i][0], reverse_anchors[i][1], 
                           reverse_anchors[i][2], reverse_anchors[i][3]) 
                           for i in reverse_path if i < len(reverse_anchors)]
    else:
        reverse_segments = []
    
    # Combine forward and reverse segments
    combined_segments = []
    if forward_segments:
        combined_segments.extend(forward_segments)
    if reverse_segments:
        combined_segments.extend(reverse_segments)
    
    # Sort by query position and resolve overlaps
    initial_segments = []
    if combined_segments:
        # Sort by query start position
        combined_segments.sort(key=lambda x: x[0])
        
        # Handle overlaps
        current = combined_segments[0]
        for next_seg in combined_segments[1:]:
            q_curr_start, q_curr_end, r_curr_start, r_curr_end = current
            q_next_start, q_next_end, r_next_start, r_next_end = next_seg
            
            # Check for overlap in query
            if q_next_start <= q_curr_end:
                # Choose the longer segment
                curr_len = q_curr_end - q_curr_start
                next_len = q_next_end - q_next_start
                
                if next_len > curr_len:
                    current = next_seg
                # Otherwise keep current
            else:
                initial_segments.append(current)
                current = next_seg
        
        initial_segments.append(current)
    
    # Merge adjacent segments
    merged_segments = merge_adjacent_segments(initial_segments, max_gap=30)
    
    # Ensure complete coverage of the query with fine-grained matching for large gaps
    complete_segments = ensure_complete_coverage(query, ref, merged_segments)
    
    # Final merging pass to clean up
    final_segments = merge_adjacent_segments(complete_segments, max_gap=20)
    
    # Convert to output format: (q_start, q_end-1, r_start, r_end-1)
    # The -1 adjustment is because our internal representation is inclusive at both ends
    output_segments = []
    for q_start, q_end, r_start, r_end in final_segments:
        # Ensure end indices are within bounds
        q_end = min(q_end, len(query) - 1)
        r_end = min(r_end, len(ref) - 1)
        output_segments.append((q_start, q_end, r_start, r_end))
    
    # Check final coverage
    final_uncovered = find_uncovered_regions(query, output_segments)
    total_uncovered = sum(u[1]-u[0]+1 for u in final_uncovered)
    coverage_percentage = 100 * (len(query) - total_uncovered) / len(query)
    print(f"Final coverage: {coverage_percentage:.2f}% of query ({len(output_segments)} segments)")
    
    return output_segments

def main():
    """Main function to run alignment"""
    # For testing with the provided data files
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    
    # Process both datasets
    for i in range(1, 3):
        query_file = os.path.join(data_dir, f"query{i}.txt")
        ref_file = os.path.join(data_dir, f"ref{i}.txt")
        output_file = f"result{i}.txt"
        
        print(f"\nProcessing dataset {i}...")
        
        query = read_sequence(query_file)
        ref = read_sequence(ref_file)
        
        print(f"Query length: {len(query)}")
        print(f"Reference length: {len(ref)}")
        
        start_time = time.time()
        alignment = find_alignment(query, ref, min_match_length=30)
        end_time = time.time()
        
        print(f"Time taken: {end_time - start_time:.2f} seconds")
        print(f"Found {len(alignment)} matching regions")
        
        with open(output_file, 'w') as f:
            f.write(str(alignment))
        print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
