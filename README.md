# Hantavirus-clockrate
This repository contains the metadata, files and scripts for running multiple analysis.

## Contents

- [Metadata](https://github.com/charu3618/Hantavirus-clockrate/tree/main/metadata): Contains curated metadata for samples used for the analysis
    - Contains metadata for S, M and L segments
    - Contains matched samples across all segments
- [ancestral reconstruction](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction): Contains files and [scripts](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/scripts) for running ancestral reconstruction using IQ-TREE
    - [Alignments](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/curated_alignments) and [rooted trees](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/rooted%20trees) were used to reconstruct ancestral sequences. Command used:
      
      ```bash
      iqtree2 -s ./sequence.fasta -te ./rooted_tree.nwk -m UNREST -asr
      ```
    - [Output files](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/iqtree_output) can be used to generate [reconstructed ancestor sequences](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/ancestral_sequences) using this [script](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/scripts)
    - These reconstructed sequences are aligned against reference sequence to get [aligned ancestral sequences](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/ancestral_sequences_aligned)
    - The [script](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/scripts) gives the final list of [mutations](https://github.com/charu3618/Hantavirus-clockrate/tree/main/ancestral%20reconstruction/mutations%20list)
- [PoW](https://github.com/charu3618/Hantavirus-clockrate/tree/main/PoW): Contains files and [scripts](https://github.com/charu3618/Hantavirus-clockrate/tree/main/PoW/scripts) for running the PoW model
    - [Rate posterior distributions](https://github.com/charu3618/Hantavirus-clockrate/tree/main/PoW/Rate_posterior_distributions): For getting posteriors of clockRate
    - [Ultrametric trees](https://github.com/charu3618/Hantavirus-clockrate/tree/main/PoW/ultrametric_trees)
    - Script output:
        - [Recreated trees](https://github.com/charu3618/Hantavirus-clockrate/tree/main/PoW/recreated_trees)
        - [Transformed trees](https://github.com/charu3618/Hantavirus-clockrate/tree/main/PoW/PoW_transformed_trees): Contains PoW transformed tree for each segment and the summarised MCC trees

- [iSNV analysis](https://github.com/charu3618/Hantavirus-clockrate/tree/main/isnv_analysis)
    -From SRA, get bam files
      PP_006XDHK.5, BioSample ID SAMN60423882;
      PP_006XDJH.5, NOT on SRA yet

    ```bash
    prefetch SRR38840885
    sam-dump SRR38840885 | samtools view -b -o SRR38840885.bam
    ```
    
    -Generate separate bam files for each segment
  
     ```bash
     samtools view -b SRR38840885.bam MT956619 > XDHK_M.bam
     samtools mpileup --no-BAQ --min-BQ 5 --max-depth 1000000 --reference ../MT956618.fasta XDHK_S.bam > XDHK_S.pileup
     ```
    -Get basefreqs

      python ~/shiver/bin/tools/AnalysePileup.py XDHK_M.pileup ../MT956619.fasta > XDHK_M_Basefreqs.csv
      
      
## Softwares used

- R v4.5.1
- BEAST v2.6.7
- IQ-TREE v2.3.5
- MAFFT v7.52
- shiver v1.7.3
- samtools v1.23.1
- python 3.14.4
  
