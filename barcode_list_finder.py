#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import random
import math


def get_seq_shingles(raw_barcodes, shingle_size):
    barcode_shingles = [set() for i in range(len(raw_barcodes))]
    all_shingles = set()
    for i, barcode in enumerate(raw_barcodes):
        for j in range(len(barcode) - shingle_size):
            shingle = barcode[j : j + shingle_size]
            all_shingles.add(shingle)
            barcode_shingles[i].add(shingle)
    return barcode_shingles, all_shingles


def hash_function(x, a, b, prime=10 ** 9 + 7):
    return (a * x + b) % prime


def get_encoded_shingle(shingle):
    encoding = 0

    for i, letter in enumerate(reversed(shingle)):
        encoding += ord(letter) << (i * 8)  # 8 bits per letter
    return encoding


def get_min_hash_matrix(
    num_of_barcodes, barcode_shingles, shingles, num_hashes, bit_hashing_depth
):
    """
    min_hash matrix is a num_of_barcodes * num_hashes matrix that stores the
    min(hash(shingle)) for all shingles in a barcode and for all the barcodes
    in a dataset. This simulates a min_valuue of permutation.
    
    """
    max_val = 2 ** bit_hashing_depth - 1
    perms = [
        [random.randint(0, max_val), random.randint(0, max_val)]
        for _ in range(num_hashes)
    ]
    min_hash_matrix = [
        [math.inf for _ in range(num_hashes)] for _ in range(num_of_barcodes)
    ]

    for i, (a, b) in enumerate(perms):
        hashes = {}
        for shingle in shingles:
            encoded_shingle = get_encoded_shingle(shingle)
            hashes[shingle] = hash_function(encoded_shingle, a, b)

        for j, barcode_set in enumerate(barcode_shingles):
            for shingle in barcode_set:
                min_hash_matrix[j][i] = min(min_hash_matrix[j][i], hashes[shingle])

    return min_hash_matrix


def get_band_code(hashes, bit_hashing_depth):
    code = 0
    for i, hashed_val in enumerate(hashes):
        code += hashed_val << (bit_hashing_depth * i)
    return code


def get_LSH_buckets(min_hash_matrix, num_hashes, bit_hashing_depth, rows_per_band):
    """
    Divide the min_hash_matrix in bands of size rows_per_band and find a hashed value
    for each barcode band. The index of the barcode is saved onto a bucket set with
    key = hashed_value. If two indexes are saved onto a same hashed_value, we have most
    likely found similarity!
    
    if bands is small, more oportunity for finding less similar sequences
    
    best t is 1/b**1/r
    """
    max_val = 2 ** bit_hashing_depth - 1
    assert (
        not num_hashes % rows_per_band
    ), "num_hashes must be divisible by rows_per_band"
    num_of_bands = num_hashes // rows_per_band
    a = [random.randint(0, max_val) for _ in range(num_of_bands)]
    b = [random.randint(0, max_val) for _ in range(num_of_bands)]

    buckets = {}

    for i, barcode_hashes in enumerate(min_hash_matrix):
        for j in range(0, num_hashes, rows_per_band):
            band_code = get_band_code(
                barcode_hashes[j : j + rows_per_band], bit_hashing_depth
            )
            bucket = hash_function(
                band_code, a[j // rows_per_band], b[j // rows_per_band]
            )
            if bucket not in buckets:
                buckets[bucket] = set()
            buckets[bucket].add(i)
    return buckets


def get_most_likely_barcode(barcode_indexes, raw_barcodes):
    """
    
    """
    barcode_len = len(raw_barcodes[0])
    nucleotides = ["A", "T", "C", "G"]
    nucleotide_count = [
        {nucleotide: 0 for nucleotide in nucleotides} for _ in range(barcode_len)
    ]

    for barcode_index in barcode_indexes:
        for i, letter in enumerate(raw_barcodes[barcode_index]):
            nucleotide_count[i][letter] += 1
    barcode = ""
    for nucleotide_i in nucleotide_count:
        max_nucleotide_count = 0
        max_nucleotide = ""
        for nucleotide, count in nucleotide_i.items():
            if count > max_nucleotide_count:
                max_nucleotide_count = count
                max_nucleotide = nucleotide
        barcode += max_nucleotide

    assert (
        len(barcode) == barcode_len
    ), "something went wrong, missing nucleotides in barcode"
    return barcode


def get_barcodes_from_buckets(buckets, raw_barcodes, sim_elements_threshold):
    barcodes = set()

    for _, barcode_indexes in buckets.items():
        if len(barcode_indexes) >= sim_elements_threshold:
            barcode = get_most_likely_barcode(barcode_indexes, raw_barcodes)
            barcodes.add(barcode)
    return barcodes


def get_barcode_list(
    raw_barcodes,
    shingle_size=7,
    num_hashes=128,
    bit_hashing_depth=32,
    LSH_rows_per_band=8,
    sim_elements_threshold=20,
):

    assert raw_barcodes, "raw_barcodes list must not be empty"
    assert math.log(num_hashes, 2).is_integer(), "num_hashes must be a power of 2"

    barcode_shingles, all_shingles = get_seq_shingles(raw_barcodes, shingle_size)
    min_hash_matrix = get_min_hash_matrix(
        len(raw_barcodes), barcode_shingles, all_shingles, num_hashes, bit_hashing_depth
    )

    buckets = get_LSH_buckets(
        min_hash_matrix, num_hashes, bit_hashing_depth, LSH_rows_per_band
    )
    barcodes = get_barcodes_from_buckets(buckets, raw_barcodes, sim_elements_threshold)

    return barcodes
