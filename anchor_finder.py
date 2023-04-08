#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 08:35:33 2023

@author: samuel
"""


import Levenshtein


def get_similarity(anchor, seq):
    return Levenshtein.distance(anchor, seq[-len(anchor) :])


def get_anchors(data, unused_nucleotides, anchor_len):
    """
    TODO: add stats on the anchor candidates
    (closest values on most frequen value, and 2nd most frequent)
    closest values on most frequent value and 3rd most frequent
    
    if anchors are too close they are assumed to be mutations of the other anchor.
    if they are infrequent, they are assumed to be mutations.
    """
    nucleotide_freq = [{"A": 0, "T": 0, "C": 0, "G": 0} for i in range(anchor_len)]

    for i, row in data.iterrows():
        for j, nucleotide in enumerate(row.seq[-anchor_len:]):
            nucleotide_freq[j][nucleotide] += 1

    anchors = ["" for i in range(4)]
    for nucleotide_count in nucleotide_freq:
        ith_nucleotide_list = list(nucleotide_count.items())
        ith_nucleotide_list.sort(key=lambda x: x[1], reverse=True)

        for i, nucleotide in enumerate(ith_nucleotide_list):
            anchors[i] += nucleotide[0]

    return anchors
