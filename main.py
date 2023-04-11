#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
from anchor_finder import get_anchor
from raw_barcode_finder import get_raw_barcodes
from barcode_list_finder import get_barcode_list
import csv


file = input("Full full file address: ")

data = pd.DataFrame(
    pd.read_csv(file, sep="\n", header=None).values.reshape(-1, 4),
    columns=["read_id", "seq", "+", "qual"],
)


barcode_len = 30
unused_nucleotides = 2
anchor_len = len(data.iloc[0].seq) - (barcode_len + unused_nucleotides)

anchor = get_anchor(data, unused_nucleotides, anchor_len)
raw_barcodes = get_raw_barcodes(anchor, data, barcode_len)
barcode_set = get_barcode_list(raw_barcodes)


with open("output.csv", "w", newline="", encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile, delimiter=",")
    for element in barcode_set:
        writer.writerow([element])
