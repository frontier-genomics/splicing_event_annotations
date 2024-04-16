#!/bin/bash
pip install -r requirements.txt
pip install -r requirements-dev.txt

#have to remove the weird chromosomes chrUn onwards
bedtools sort -g chromosome_names.txt -i refseq_curated_introns.bed > refseq_curated_introns_sorted.bed 
bedtools sort -g chromosome_names.txt -i refseq_mane_introns.bed > refseq_mane_introns_sorted.bed