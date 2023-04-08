#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 08:38:57 2023

@author: samuel
"""


def is_valid_seq(anchors, seq, sequence_anchor_similarity_threshold=0.9):
    """ 
    TODO:
    - improve performance and add shift in anchor support 
    - if first nucleotides are not present, don't use sequence
    
    Valid sequences have a valid anchor
    """
    max_similarity = 0
    best_anchor_i = -1
    shifts = -1

    for i, anchor in enumerate(anchors):
        anchor_len = len(anchor)
        for shift in range(3):
            similarity = 0
            for (anchor_letter, seq_letter) in zip(
                anchor[shift:], seq[-anchor_len + shift :]
            ):
                if anchor_letter == seq_letter:
                    similarity += 1
            if similarity > max_similarity:
                max_similarity = similarity
                shifts = shift
                best_anchor_i = i

    if max_similarity > int(
        len(anchors[best_anchor_i]) * sequence_anchor_similarity_threshold
    ):
        return True, shifts
    else:
        return False, -1


def get_barcode(anchors, seq, barcode_len):
    is_valid, shifts = is_valid_seq(anchors, seq)
    if is_valid:
        return seq[shifts : shifts + barcode_len]
    else:
        return ""


def get_raw_barcodes(anchors, data, barcode_len):
    barcodes = []
    for _, row in data.iterrows():
        seq = row.seq
        barcode = get_barcode(anchors, seq, barcode_len)
        if barcode:
            barcodes.append(barcode)
    return barcodes
