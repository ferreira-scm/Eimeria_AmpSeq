#!/bin/bash

qiime rescript get-silva-data \
      --p-version '138.1' \
      --p-target 'LSURef_NR99' \
      --p-include-species-labels \
      --p-ranks kingdom phylum class order family genus \
      --o-silva-sequences silva-138.1-lsu-nr99-rna-seqs.qza \
      --o-silva-taxonomy silva-138.1-lsu-nr99-tax.qza

qiime rescript reverse-transcribe \
      --i-rna-sequences silva-138.1-lsu-nr99-rna-seqs.qza \
      --o-dna-sequences silva-138.1-lsu-nr99-seqs.qza

qiime rescript cull-seqs \
      --i-sequences silva-138.1-lsu-nr99-seqs.qza \
      --o-clean-sequences silva-138.1-lsu-nr99-seqs-cleaned.qza

#qiime rescript filter-seqs-length-by-taxon \
#      --i-sequences  silva-138.1-lsu-nr99-rna-seqs-cleaned.qza\
#      --i-taxonomy silva-138.1-lsu-nr99-tax.qza \
#      --p-labels Archaea Bacteria Eukaryota \
#      --p-min-lens 2000 2000 4000 \
#      --o-filtered-seqs silva-138.1-lsu-nr99-rna-seqs-fil.qza \
#      --o-discarded-seqs silva-138.1-lsu-nr99-seqs-discard.qza 

qiime rescript dereplicate \
      --i-sequences silva-138.1-lsu-nr99-seqs-cleaned.qza \
      --i-taxa silva-138.1-lsu-nr99-tax.qza \
      --p-rank-handles 'silva' \
      --p-mode 'majority' \
      --o-dereplicated-sequences silva-138.1-lsu-nr99-seqs-derep-maj.qza \
      --o-dereplicated-taxa silva-138.1-lsu-nr99-tax-derep-maj.qza

qiime tools export \
      --input-path silva-138.1-lsu-nr99-tax-derep-maj.qza \
      --output-path Fastas/

qiime tools export \
      --input-path silva-138.1-lsu-nr99-seqs-derep-maj.qza \
      --output-path Fastas/

mv Fastas/dna-sequences.fasta Fastas/silva-138.1-lsu-nr99-seqs.fasta
mv Fastas/taxonomy.tsv Fastas/silva-138.1-lsu-nr99-tax.tsv

qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads silva-138.1-lsu-nr99-seqs-derep-maj.qza \
      --i-reference-taxonomy silva-138.1-lsu-nr99-tax-derep-maj.qza \
      --o-classifier silva-138.1-lsu-nr99-classifier.qza

qiime tools export \
      --input-path silva-138.1-lsu-nr99-classifier.qza \
      --output-path Fastas/
