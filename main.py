#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:14:17 2023

@author: samuel

Assumption:
    sequences only include "A" "T" "G" or "C" nucleotides

Complexity:len
    num_of_raw_barcodes * len_barcode * num_hash_functions  + num_of_extracted

Limitations: only supports 1 anchor, however it could be improved using
            the barcode retrieval method 
"""

import pandas as pd
from anchor_finder import get_anchor
from raw_barcode_finder import get_raw_barcodes
from barcode_list_finder import get_barcode_list


data = pd.DataFrame(
    pd.read_csv(
        "/home/samuel/Documents/Work/Cajal/tail_100000.fastq", sep="\n", header=None
    ).values.reshape(-1, 4),
    columns=["read_id", "seq", "+", "qual"],
)


barcode_len = 30
unused_nucleotides = 2
anchor_len = len(data.iloc[0].seq) - (barcode_len + unused_nucleotides)

anchor = get_anchor(data, unused_nucleotides, anchor_len)
raw_barcodes = get_raw_barcodes(anchor, data, barcode_len)
barcode_set = get_barcode_list(raw_barcodes)
