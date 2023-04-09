#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:14:17 2023

@author: samuel

Assumption:
    sequences only include "A" "T" "G" or "C" nucleotides

"""

import pandas as pd
from anchor_finder import get_anchors
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

anchors = get_anchors(data, unused_nucleotides, anchor_len)
raw_barcodes = get_raw_barcodes(anchors, data, barcode_len)
barcode_list = get_barcode_list(raw_barcodes)
