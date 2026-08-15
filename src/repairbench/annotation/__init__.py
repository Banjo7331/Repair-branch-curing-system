"""Reading real reference annotation instead of assuming it.

Everything above this package has been reasoning about transcripts whose
structure was asserted in a fixture. That was fine for testing rules and is not
fine for anything else: the NMD calculation — which decides between two opposite
therapies — is arithmetic over exon boundaries, and arithmetic over made-up
boundaries is decoration.

This package turns asserted inputs into earned ones. It parses GFF3, reads a
reference FASTA, maps genomic coordinates onto coding positions, and left-aligns
indels the way the rest of the field does. Each source is pinned by content
digest, so a call can name the annotation release it was made against.
"""
