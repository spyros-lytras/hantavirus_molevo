# Hantavirus-clockrate
This repository contains the metadata, files and scripts for running multiple analysis for *Molecular evolution of Andes virus lineages leading to recent human outbreaks*.

## Contents

- [Metadata](metadata): Contains curated metadata for samples used for the analysis
    - Contains metadata for S, M and L segments
    - Contains matched samples across all segments
- [ancestral reconstruction](ancestral_reconstruction): Contains files and [scripts](ancestral_reconstruction/scripts) for running ancestral reconstruction using IQ-TREE
    - [Alignments](codon_alignments) and [rooted trees](phylogenies) were used to reconstruct ancestral sequences. Command used:
      
      ```bash
      iqtree2 -s ./sequence.fasta -te ./rooted_tree.nwk -m UNREST -asr
      ```
    - [Output files](ancestral_reconstruction/iqtree_output) can be used to generate [reconstructed ancestor sequences](ancestral_reconstruction/ancestral_sequences) using this [script](ancestral_reconstruction/scripts)
    - These reconstructed sequences are aligned against reference sequence to get [aligned ancestral sequences](ancestral_reconstruction/ancestral_sequences_aligned)
    - The [script](ancestral_reconstruction/scripts) gives the final list of [mutations](ancestral_reconstruction/mutations_list)
- [Beast analysis](beast): Contains BEAST2 XML files and the random number seeds used to run the molecular dating analyses.
- [BETS] (https://github.com/spyros-lytras/hantavirus_molevo/tree/main/BETS): Contains BEAST XMLs for running BETS analysis on strict clock and relaxed clock datasets, including/excluding starting tree.
- [PoW](PoW): Contains files and [scripts](PoW/scripts) for running the PoW model
    - [Rate posterior distributions](PoW/Rate_posterior_distributions): For getting posteriors of clockRate
    - [Ultrametric trees](PoW/ultrametric_trees)
    - Script output:
        - [Recreated trees](PoW/recreated_trees)
        - [Transformed trees](PoW_transformed_trees): Contains PoW transformed tree for each segment and the summarised MCC trees

- [iSNV analysis](isnv_analysis)
    -From SRA, get bam files
      PP_006XDHK.5, BioSample ID SAMN60423882;
      PP_006XDJH.5, BioSample ID SAMN60696784

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

- [Selection analysis](selection_analysis): Contains all scripts, documentation and results for the selection analysis.    
      
## Softwares used

- R v4.5.1
- BEAST v2.6.7
- IQ-TREE v2.3.5
- MAFFT v7.52
- shiver v1.7.3
- samtools v1.23.1
- python 3.14.4
  
