#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 08:49:02 2023

@author: samuel
"""

import numpy as np
import random
import pandas as pd


def get_random_seq(seq_len):
    nucleotides = ["A", "T", "C", "G"]
    seq = ""
    for _ in range(seq_len):
        seq += random.choice(nucleotides)
    return seq


def mutate_seq(current_seq, mutation_type, nucleotides):

    mutation_index = random.randint(0, len(current_seq) - 1)
    if mutation_type == "add":
        current_seq = (
            current_seq[:mutation_index]
            + random.choice(nucleotides)
            + current_seq[mutation_index:]
        )
    elif mutation_type == "delete":
        current_seq = current_seq[:mutation_index] + current_seq[mutation_index + 1 :]
    elif mutation_type == "subsitute":
        current_seq = (
            current_seq[:mutation_index]
            + random.choice(nucleotides)
            + current_seq[mutation_index + 1 :]
        )
    return current_seq


def get_test_sequences(
    num_of_anchors,
    anchor_probability,
    barcode_len,
    unused_nucleotides,
    anchor_len,
    num_of_sequences,
    mutation_probability,
    num_of_barcodes,
    mutation_type,
):
    """
    Generate random sequences to test if anchor_finder finds the correct anchor.
    
    This function uses a random mutation probability to shift the anchor and to
    simulate errors in transcription.
    """

    anchors = [get_random_seq(anchor_len) for _ in range(num_of_anchors)]
    barcodes = [get_random_seq(barcode_len) for _ in range(num_of_barcodes)]

    nucleotides = ["A", "T", "C", "G"]
    sequence_size = barcode_len + unused_nucleotides + anchor_len

    sequences = []
    for _ in range(num_of_sequences):
        current_anchor = np.random.choice(anchors, 1, p=anchor_probability)[0]

        # Mutate
        while (
            random.random() < mutation_probability
            and 0 < len(current_anchor) < sequence_size
        ):
            current_anchor = mutate_seq(current_anchor, mutation_type, nucleotides)

        current_barcode = random.choice(barcodes)

        while (
            random.random() < mutation_probability
            and 0 < len(current_barcode) + len(current_anchor) < sequence_size
        ):
            current_barcode = mutate_seq(current_barcode, mutation_type, nucleotides)

        current_unused_nucleotides = "".join(
            [random.choice(nucleotides) for _ in range(2)]
        )
        current_sequence = current_barcode + current_unused_nucleotides + current_anchor

        while len(current_sequence) > sequence_size:
            del_nucleotide_i = random.randint(0, len(current_sequence) - 1)
            current_sequence = (
                current_sequence[:del_nucleotide_i]
                + current_sequence[del_nucleotide_i + 1 :]
            )
        while len(current_sequence) < sequence_size:
            del_nucleotide_i = random.randint(0, len(current_sequence) - 1)
            current_sequence = (
                current_sequence[:del_nucleotide_i]
                + random.choice(nucleotides)
                + current_sequence[del_nucleotide_i:]
            )

        assert (
            len(current_sequence) == sequence_size
        ), "incorrect sequence size produced"
        sequences.append(current_sequence)
    return pd.DataFrame(sequences, columns=["seq"]), anchors, barcodes
