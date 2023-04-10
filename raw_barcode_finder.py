#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 08:38:57 2023

@author: samuel
"""


def get_shingles(seq, shingle_size=3):
    assert len(seq) > shingle_size, "sequence is too short, decrease shingle_size"
    shingles = set()

    for i in range(0, len(seq) - shingle_size):
        shingles.add(seq[i : i + shingle_size])
    return shingles


def jaccard_similarity(seq, anchor):
    anchor_len = len(anchor)
    seq_shingles = get_shingles(seq[-anchor_len:])
    anchor_shingles = get_shingles(anchor)

    return len(seq_shingles.intersection(anchor_shingles)) / len(
        seq_shingles.union(anchor_shingles)
    )


def get_shifts(anchor, seq, minimum_identical_start):
    anchor_start = anchor[:minimum_identical_start]
    barcode_len = len(seq) - len(anchor)
    assert barcode_len > 0, "something went wrong anchor retrieved is too long"

    for shifts in range(len(anchor) - minimum_identical_start):
        start = shifts + barcode_len
        end = shifts + barcode_len + minimum_identical_start
        if seq[start:end] == anchor_start:
            return shifts
    return -1


def is_valid_seq(
    anchor, seq, seq_anchor_similarity_threshold=0.7, min_identical_start=5
):
    """
    For a seq to have a valid barcode it must:
        1. Have an anchor that starts at  index barcode_len+unused_nucleotides_len+1 so there is no loss
           of barcode nucleotides. 
        2. Have an anchor at least "seq_anchor_similarity_threshold" similar to the actual anchor.
        3. Start with at least "min_identical_start" nucleotides. This improves
           the barcode quality by avoiding mutations where nucleotides
           could have been added, deleted or substituted, which would retrieve
           different shifted versions of the barcode. 
           3.1 This rule could be changed by obtaining the best jaccard_similarity
               between seq_anchor shifted versions and anchor and accepting
               possible errors in the barcode. This is slightly slower.
    """

    seq_anchor_similarity = jaccard_similarity(seq, anchor)
    if seq_anchor_similarity >= seq_anchor_similarity_threshold:
        shifts = get_shifts(anchor, seq, min_identical_start)
        if shifts != -1:
            return True, shifts
    return False, -1


def get_barcode(anchor, seq, barcode_len):
    is_valid, shifts = is_valid_seq(anchor, seq)
    if is_valid:
        return seq[shifts : shifts + barcode_len]
    else:
        return ""


def get_raw_barcodes(anchor, data, barcode_len):
    barcodes = []
    for _, row in data.iterrows():
        seq = row.seq
        barcode = get_barcode(anchor, seq, barcode_len)
        if barcode:
            barcodes.append(barcode)
    return barcodes
