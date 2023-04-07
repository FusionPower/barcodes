#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:59:17 2023

@author: samuel
"""

from anchor_finder import get_anchors
import pandas as pd
import random


def get_random_anchor(anchor_len):
    nucleotides = ["A", "T", "C", "G"]
    anchor = ""
    for _ in range(anchor_len):
        anchor += random.choice(nucleotides)
    return anchor


def get_test_sequences(barcode_len, unused_nucleotides, anchor_len, num_of_sequences, mutation_probability, anchor_mutation_type):
    """
    Generate random sequences to test if anchor_finder finds the correct anchor.
    
    This function uses a random mutation probability to shift the anchor and to
    simulate errors in transcription.
    """
    nucleotides = ["A", "T", "C", "G"]
    sequence_size = barcode_len + unused_nucleotides + anchor_len
    anchor = get_random_anchor(anchor_len)

    sequences = []
    for _ in range(num_of_sequences):
        current_sequence = anchor
        
        # Mutate
        while random.random() < mutation_probability and 0 < len(current_sequence) < sequence_size:
            shift_index = random.randint(0, anchor_len-1)
            if anchor_mutation_type == "shift":
                current_sequence = anchor[:shift_index] + random.choice(nucleotides) + anchor[shift_index:]
            elif anchor_mutation_type == "delete":
                current_sequence = anchor[:shift_index] + anchor[shift_index+1:]
            elif anchor_mutation_type == "subsitute":
                current_sequence = anchor[:shift_index] + random.choice(nucleotides) + anchor[shift_index+1:]
        
        for _ in range(sequence_size - len(current_sequence)):
            current_sequence = random.choice(nucleotides) + current_sequence
        sequences.append(current_sequence)
    return pd.DataFrame(sequences, columns=["seq"]), anchor



def test_anchor_finder(
    barcode_len=30,
    unused_nucleotides=2,
    anchor_len=47,
    num_of_sequences=10000,
    mutation_anchor_probability=0.2,
):

    sequences, real_anchor = get_test_sequences(
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_anchor_probability,
        anchor_mutation_type = "shift"
    )
    found_anchors = get_anchors(sequences, unused_nucleotides, anchor_len)
    
    assert real_anchor == found_anchors[0]

    sequences, real_anchor = get_test_sequences(
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_anchor_probability,
        anchor_mutation_type = "subtitute"
    )
    found_anchors = get_anchors(sequences, unused_nucleotides, anchor_len)
    
    assert real_anchor == found_anchors[0]
    
    sequences, real_anchor = get_test_sequences(
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_anchor_probability,
        anchor_mutation_type = "delete"
    )
    found_anchors = get_anchors(sequences, unused_nucleotides, anchor_len)
    
    assert real_anchor == found_anchors[0]