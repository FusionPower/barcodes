#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 05:18:11 2023

@author: samuel
"""

import random
import math


def get_seq_shingles(raw_barcodes, shingle_size):
    barcode_shingles = [set() for i in range(len(raw_barcodes))]
    shingles = set()
    for i, barcode in enumerate(raw_barcodes):
        for j in range(len(barcode) - shingle_size):
            shingle = barcode[j : j + shingle_size + 1]
            shingles.add(shingle)
            barcode_shingles[i].add(shingle)
    return barcode_shingles, shingles


def hash_function(x, a, b, prime=10 ** 9 + 7):
    return (a * x + b) % prime


def get_encoded_shingle(shingle):
    encoding = 0
    for i, letter in enumerate(reversed(shingle)):
        encoding += ord(letter) << i
    return encoding


def get_min_hash_matrix(
    len_barcodes, barcode_shingles, shingles, num_hashes, bit_hashing_depth
):
    max_val = 2 ** bit_hashing_depth - 1
    perms = [
        [random.randint(0, max_val), random.randint(0, max_val)]
        for _ in range(num_hashes)
    ]
    min_hash_matrix = [
        [math.inf for _ in range(num_hashes)] for _ in range(len_barcodes)
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


def get_LSH_buckets(min_hash_matrix, num_hashes, bit_hashing_depth):
    """
    TODO: move r to a more accesible function
    """
    max_val = 2 ** bit_hashing_depth - 1
    a, b = random.randint(0, max_val), random.randint(0, max_val)
    r = 8

    buckets = {}

    for i,barcode_hashes in enumerate(min_hash_matrix):
        for j in range(0, num_hashes, r):
            band_code = get_band_code(barcode_hashes[j : j + r], bit_hashing_depth)
            bucket = hash_function(band_code, a, b)
            if bucket not in buckets:
                buckets[bucket] = set()
            buckets[bucket].add(i)
    return buckets


def get_most_likely_barcode(barcode_indexes, raw_barcodes):
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


def get_barcodes_from_buckets(buckets, raw_barcodes):
    barcodes = set()

    for _, barcode_indexes in buckets.items():
        if len(barcode_indexes) > 1:
            barcode = get_most_likely_barcode(barcode_indexes, raw_barcodes)
            barcodes.add(barcode)
    return barcodes


def get_barcode_list(
    raw_barcodes, shingle_size=7, num_hashes=128, bit_hashing_depth=32
):

    assert raw_barcodes, "raw_barcodes list must not be empty"
    assert math.log(num_hashes, 2).is_integer(), "num_hashes must be a power of 2"
    
    barcode_shingles, shingles = get_seq_shingles(raw_barcodes, shingle_size)
    min_hash_matrix = get_min_hash_matrix(
        len(raw_barcodes), barcode_shingles, shingles, num_hashes, bit_hashing_depth
    )

    buckets = get_LSH_buckets(min_hash_matrix, num_hashes, bit_hashing_depth)
    barcodes = get_barcodes_from_buckets(buckets, raw_barcodes)

    return barcodes
