### Scripts here are used to preprocess sequences and extract indels from genomic sequences.

### Usage:

1. preprocess: Install nextstrain environment according to https://docs.nextstrain.org/en/latest/install.html

```bash
conda activate nextstrain

# Sanitize metadata.
python3 sanitize_metadata.py --metadata metadata_tsv_2022_08_18.tar.xz --database-id-columns "Accession ID" --parse-location-field Location --rename-fields 'Virus name=strain' 'Accession ID=gisaid_epi_isl' 'Collection date=date' --strip-prefixes "hCoV-19/" --output 20220818_metadata.rmdup.tsv.gz

# Sanitize sequences
python3 sanitize_sequences.py --sequences sequences_fasta_2022_08_18.tar.xz --strip-prefixes "hCoV-19/" --output 20220818_sequences.rmdup.fa.gz

# Generate SAM file
minimap2 -a -x asm5 -t 36 NC_045512.2.fasta 20220818_sequences.rmdup.fa.gz > 01_sequences.ref.sam
```

2. extract indels with scripts **indel_identification.py**:

```bash
python indel_identification.py --sam 01_sequences.ref.sam --output 02_indel.info.tsv --summary 03_num_indel.seq.tsv
```
