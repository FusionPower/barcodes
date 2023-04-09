#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:59:17 2023

@author: samuel
"""

from anchor_finder import get_anchors
from testing_utilities import get_test_sequences


def test_anchor_finder(
    barcode_len=30,
    unused_nucleotides=2,
    anchor_len=47,
    num_of_sequences=50000,
    mutation_anchor_probability=0.2,
    num_of_anchors=4,
):
    anchor_probability = [0.8, 0.1, 0.05, 0.05]
    assert (
        len(anchor_probability) == num_of_anchors
    ), "probability per anchor must match number of anchors"

    # Test shift mutations
    sequences, real_anchors, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_anchor_probability,
        mutation_type="add",
    )

    found_anchors = get_anchors(sequences, unused_nucleotides, anchor_len)
    assert (
        real_anchors[0] == found_anchors[0]
    ), "real anchor not found for mutation 'add'"

    # Test substitute mutations
    sequences, real_anchors, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_anchor_probability,
        mutation_type="subtitute",
    )
    found_anchors = get_anchors(sequences, unused_nucleotides, anchor_len)
    assert (
        real_anchors[0] == found_anchors[0]
    ), "real anchor not found for mutation 'substitute'"

    # Teste delete mutations
    sequences, real_anchors, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_anchor_probability,
        mutation_type="delete",
    )
    found_anchors = get_anchors(sequences, unused_nucleotides, anchor_len)
    assert (
        real_anchors[0] == found_anchors[0]
    ), "real anchor not found for mutation 'delete'"
