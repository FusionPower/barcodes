#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 08:46:19 2023

@author: samuel
"""


from barcode_list_finder import get_barcode_list
from raw_barcode_finder import get_raw_barcodes
from testing_utilities import get_test_sequences


def test_barcode_finder(
    barcode_len=30,
    unused_nucleotides=2,
    anchor_len=47,
    num_of_sequences=10,
    mutation_probability=0.5,
    num_of_anchors=1,
    num_of_barcodes=1,
):

    anchor_probability = [1]  # Probability distribution for anchors
    assert (
        len(anchor_probability) == num_of_anchors
    ), "probability per anchor must match number of anchors"
    assert num_of_anchors == 1, "num_of_anchors only supports size one anchor for how"

    # TEST ADD MUTATIONS
    sequences, anchor, barcodes = get_test_sequences(
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

    raw_barcodes = get_raw_barcodes(anchor, sequences, barcode_len)
    barcode_set = get_barcode_list(raw_barcodes)

    jaccard_similarity = len(barcode_set) / len(barcode_set.union(barcodes))
    print(f"Jaccard_similarity = {jaccard_similarity}")
    assert jaccard_similarity > 0.95, "Original barcodes were not found"

    # TEST DELETE MUTATIONS
    sequences, anchor, barcodes = get_test_sequences(
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

    raw_barcodes = get_raw_barcodes(anchor, sequences, barcode_len)
    barcode_set = get_barcode_list(raw_barcodes)

    jaccard_similarity = len(barcode_set) / len(barcode_set.union(barcodes))
    print(f"Jaccard_similarity = {jaccard_similarity}")
    assert jaccard_similarity > 0.95, "original barcodes were not found"

    # TEST SUBSTITUTE MUTATIONS
    sequences, anchor, barcodes = get_test_sequences(
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

    raw_barcodes = get_raw_barcodes(anchor, sequences, barcode_len)
    barcode_set = get_barcode_list(raw_barcodes)

    jaccard_similarity = len(barcode_set) / len(barcode_set.union(barcodes))
    print(f"Jaccard_similarity = {jaccard_similarity}")
    assert jaccard_similarity > 0.95, "original barcodes were not found"

    # TEST ANY MUTATION TYPE
    sequences, anchor, barcodes = get_test_sequences(
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

    raw_barcodes = get_raw_barcodes(anchor, sequences, barcode_len)
    barcode_set = get_barcode_list(raw_barcodes)

    jaccard_similarity = len(barcode_set) / len(barcode_set.union(barcodes))
    print(f"Jaccard_similarity = {jaccard_similarity}")
    assert jaccard_similarity > 0.95, "original barcodes were not found"
