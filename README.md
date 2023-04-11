# Barcodes

[![Python application](https://github.com/FusionPower/barcodes/actions/workflows/python-app.yml/badge.svg)](https://github.com/FusionPower/barcodes/actions/workflows/python-app.yml)


## Algorithm Description
The main idea used to solve this problem is Local Sensitive Hashing with Min-
hashing. The following explains the details of the implementation.
### Main module
The logic behind the main module is very straightforward:
• Extract the sequences from the data
• Retrieve the anchor with the”get anchor’ function
• Use the retrieved anchor to get the raw barcodes from the sequences us-
ing”get raw barcodes function’
• Use retrieved raw barcodes to get the shortened list of barcodes.
### Data Extraction
The extraction was done with pandas and the data was assumed to be uncor-
rupted. This means there is an iterative structure that follows the following
order ”read id”, ”sequence”, ”+”, ”quality”. Deviations from this structure
will cause faulty output or crashes.
Module Time Complexity O(number of sequences * sequence length)
## Anchor Finder module
This module strives to take advantage of the fact that, most of the time, a
single anchor is used to extract sequences. The module counts the most common
nucleotide per anchor position. It is expected that a random mutation would
yield anomalous sequences. However, the most common nucleotide per position
would be a part of the original anchor transcription. The barcode length and
the unused nucleotides between the barcode and the anchor give the first anchor
position. It is worth noting that the most detrimental mutations for this module
are shifts in the anchor since it means a lot less indexes of the true anchor will
match with the anchor candidates. A justification for the validity of this strategy
is explored in a later section.
This module only finds one anchor; however it could be improved to get-
ting up to four anchors in time complexity O(n*m) where n is the number of
sequences and m is the anchor length. This improvement is only possible if the
probability of each anchor is not uniform. As the probability of each anchor
appearing gets closer, this idea fails to retrieve a valid anchor.
If more than four anchors were needed the module logic would have to be
changed and the same idea used to shorten the list of barcodes could be used
to extract the anchors.
Another improvement possibility could be to take into account the quality
of the sequence and discard the nucleotides that have poorer quality.
It is worth mentioning that if the anchor length is known for this data, this
module could be bypassed, and the first barcode elements could be retrieved
directly. This would only slightly decrease the raw barcodes’ quality but increase
the performance.
Module Time Complexity: O(anchor length * number of sequences)

### Raw Barcode Finder module
To extract the raw barcodes from the anchor, the similarity between the se-
quence at the anchor positions and the true anchor are obtained using the
Jaccard Similarity, which is defined as:

$$ |S 1 ∩ S 2 | |S 1 ∪ S 2 | $$



Each set is composed of shingles from the original sequence. Shingles are
obtained by extracting all possible continuous substrings of size K from the
original sequence at the anchor positions. Sets are compared, and a similarity
threshold is used to discriminate between valid and non-valid anchor candidates.
This way, anchor candidates that are similar to the real anchor are interpreted
as valid anchors.
After an anchor is found to be valid, the next thing to know is to find whether
there are any shifts on the anchor such that the barcode starts at”shift’ positions
from the 0th index. To do this, there must be an exact match between the real
anchor’s first E elements and the anchor candidate’s first E elements. This
is revised starting at the first position of the possible anchor start and then
moved onward until an exact match is found. If the beginning of the anchor
is corrupted, it is impossible to know if there was a shift that would make
the raw barcode incomplete. There could be a mutation in the unused middle
nucleotides, but there is no way to tell if those were mutated with additions or
deletions since they don’t have a predefined value.
Module Time Complexity: O(sequence length * number of sequences)

