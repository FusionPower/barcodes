#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:59:17 2023

@author: samuel
"""

from anchor_finder import get_anchor
from testing_utilities import get_test_sequences


def test_anchor_finder(
    barcode_len=30,
    unused_nucleotides=2,
    anchor_len=47,
    num_of_sequences=50000,
    mutation_probability=0.5,
    num_of_anchors=1,
    num_of_barcodes=100,
):
    
    """
    This test looks for EXACT matches between real_anchor and found_anchor, however, in practice
    a similar anchor is good enough and exact matches are only important at the start of the anchor
    """

    anchor_probability = [1]
    assert (
        len(anchor_probability) == num_of_anchors
    ), "probability per anchor must match number of anchors"
    assert num_of_anchors == 1, "only anchor_len is suported for now"

    # Test shift mutations
    sequences, real_anchor, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_probability,
        num_of_barcodes,
        mutation_type="add",
    )

    found_anchor = get_anchor(sequences, unused_nucleotides, anchor_len)
    assert real_anchor == found_anchor, "real anchor not found for mutation 'add'"

    # Test substitute mutations
    sequences, real_anchor, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_probability,
        num_of_barcodes,
        mutation_type="substitute",
    )
    found_anchor = get_anchor(sequences, unused_nucleotides, anchor_len)
    assert (
        real_anchor == found_anchor
    ), "real anchor not found for mutation 'substitute'"

    # Teste delete mutations
    sequences, real_anchor, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_probability,
        num_of_barcodes,
        mutation_type="delete",
    )
    found_anchor = get_anchor(sequences, unused_nucleotides, anchor_len)
    assert real_anchor == found_anchor, "real anchor not found for mutation 'delete'"

    # Teste delete mutations
    sequences, real_anchor, _ = get_test_sequences(
        num_of_anchors,
        anchor_probability,
        barcode_len,
        unused_nucleotides,
        anchor_len,
        num_of_sequences,
        mutation_probability,
        num_of_barcodes,
        mutation_type="any",
    )
    found_anchor = get_anchor(sequences, unused_nucleotides, anchor_len)
    assert (
        real_anchor == found_anchor
    ), "real anchor not found for random mutation of type 'any'"
