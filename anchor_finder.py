#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def get_anchor(data, unused_nucleotides, anchor_len):
    """
    LOGIC:
    Gets anchor by counting the most frequent nucleotide per position and concatenating them
    
    POSSIBLE IMPROVEMENTS:
    This could be improved to getting up to 4 anchors in O(num_seq) * anchor_len:
        - first get a rough estimate of the anchors by getting the most frequent nucleotides per position and grouping them by order
        - group the sequences according to their jaccard similarity to each anchor candidate between shigles_sets of shingle size X
        - recount the most common nucleotide for each group
    
    LIMITATIONS:
        As multiple anchor probability distribution gets uniform, this algorithm crashes!
        To solve this we could use the same algorithm used to extract barcodes but performance decreases
    """

    nucleotide_freq = [{"A": 0, "T": 0, "C": 0, "G": 0} for i in range(anchor_len)]

    for _, row in data.iterrows():
        for j, nucleotide in enumerate(row.seq[-anchor_len:]):
            nucleotide_freq[j][nucleotide] += 1

    anchor = ""
    for nucleotide_count in nucleotide_freq:
        most_freq_nucleotide = ""
        max_freq = 0
        for nucleotide, count in nucleotide_count.items():
            if count > max_freq:
                max_freq = count
                most_freq_nucleotide = nucleotide
        anchor += most_freq_nucleotide
    return anchor
