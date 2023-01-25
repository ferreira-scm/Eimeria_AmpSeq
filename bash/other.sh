#!/bin/bash

qiime rescript get-ncbi-data \
      --p-n-jobs 5 \
      --p-ranks 'kingdom' 'phylum'  'class'   'order'   'family'  'genus' 'species' \
      --verbose \
      --p-query '(cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR COX1[Title] OR CO1[Title] OR COI[Title] OR MT-RNR1[Title] OR MTRNR1 12S ribosomal RNA[Title] OR 12S ribosomal RNA[Title] OR cytochrome c oxidase subunit 3[Title] OR cytochrome c oxidase subunit III[Title] OR cytochrome oxidase subunit 3[Title] OR cytochrome oxidase subunit III[Title] OR COX3[Title] OR CO3[Title] OR COIII[Title] OR mitochondrion[Title] OR Glutamate dehydrogenase[Title] OR Gdh[Title]  OR MTRNR2 16S ribosomal RNA[Title] OR MT-RNR2[Title] OR rbcL[Title] OR ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit[Title] OR  TPI[Title] OR triosephosphate isomerase[Title] OR 12S[Title] OR 12S rRNA[Title] OR (16S rRNA AND Metazoa))' \
      --output-dir NCBI_Other

qiime rescript cull-seqs \
      --i-sequences ./NCBI_Other/sequences.qza \
      --p-num-degenerates 5 \
      --o-clean-sequences NCBIdata_filtd_seqs.qza

qiime rescript filter-seqs-length \
      --i-sequences NCBIdata_filtd_seqs.qza \
      --p-global-min 300 \
      --o-filtered-seqs NCBIdata_clean_filtd_seqs.qza \
      --o-discarded-seqs NCBIdata_length_discarded_seqs.qza

qiime rescript dereplicate --verbose \
      --i-sequences NCBIdata_clean_filtd_seqs.qza \
      --i-taxa NCBI_Other/taxonomy.qza \
      --p-mode 'super' \
      --p-rank-handles 'silva' \
      --o-dereplicated-sequences NCBIdata_derep_seqs.qza \
      --o-dereplicated-taxa NCBIdata_derep_taxa.qza

qiime tools export \
      --input-path NCBIdata_derep_seqs.qza \
      --output-path Fastas/

qiime tools export \
      --input-path NCBIdata_derep_taxa.qza \
      --output-path Fastas/

mv Fastas/dna-sequences.fasta Fastas/NCBI_seqs.fasta
mv Fastas/taxonomy.tsv Fastas/NCBI_tax.tsv

#qiime feature-classifier fit-classifier-naive-bayes \
#      --i-reference-reads NCBIdata_derep_seqs.qza \
#      --i-reference-taxonomy NCBIdata_derep_taxa.qza \
#      --o-classifier NCBIdata-classifier.qza

#qiime tools export \
#      --input-path NCBIdata-classifier.qza \
#      --output-path Fastas/
