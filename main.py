#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:14:17 2023

@author: samuel
"""

import pandas as pd
from anchor_finder import get_anchors
from raw_barcode_finder import get_raw_barcodes


data = pd.DataFrame(
    pd.read_csv(
        "/home/samuel/Documents/Work/Cajal/tail_100000.fastq", sep="\n", header=None
    ).values.reshape(-1, 4),
    columns=["read_id", "seq", "+", "qual"],
)


seq_len = 30
unused_nucleotides = 2
anchor_len = len(data.iloc[0].seq) - (seq_len + unused_nucleotides)

anchors = get_anchors(data, unused_nucleotides, anchor_len)

raw_barcodes = get_raw_barcodes(data, anchors)

# shortened_barcode_list = get_shortened_barcodes
