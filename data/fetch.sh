#!/bin/sh
# Target transcripts (NCBI nuccore) and non-target transcriptomes (RefSeq RNA).
set -e
cd "$(dirname "$0")"
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
B=https://ftp.ncbi.nlm.nih.gov/genomes/refseq
curl -sS "$E/efetch.fcgi?db=nuccore&id=HQ634702,XM_049999838,XM_049985731,XM_049990443&rettype=fasta&retmode=text" -o targets.fa
curl -sS -o bee.fna.gz     "$B/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_rna.fna.gz"
curl -sS -o monarch.fna.gz "$B/invertebrate/Danaus_plexippus/latest_assembly_versions/GCF_018135715.1_MEX_DaPlex/GCF_018135715.1_MEX_DaPlex_rna.fna.gz"
curl -sS -o human.fna.gz   "$B/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna.fna.gz"
