# Computational pipeline to detect gene loss events

In order to achieve high specificity, this pipeline implements a number of steps that overcome assembly and alignment issues, and address evolutionary exon-intron structure changes in genes that are conserved.

The starting point of the pipeline is a multiple genome alignment. The output is a list of lost genes.
