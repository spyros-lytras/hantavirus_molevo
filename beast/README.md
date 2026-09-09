# BEAST2 XML FILES

This directory contains the BEAST2 XML files used for the analyses of the
S, M, and L genomic segments.

The directory is organised into three subfolders, corresponding to the
three segments:

    S_segment/
    M_segment/
    L_segment/

Each segment-specific folder contains all BEAST2 XML files used for the
respective analyses.

The analyses were performed using BEAST2 v2.6.7.

## FILE NAMING

Files follow the naming convention:

    segment_dataset_scenario.xml

where:

    segment   = genomic segment (S, M, or L)

    dataset   = clade composition:
                Clade1and3
                Clade2
                Clade1-3
                Clade1-4
                inclUnclassified

    scenario  = analysis scenario:
                allObserved
                TwoPerOutbreak
                linkageInformed

For example:

    S_Clade1-3_linkageInformed.xml

corresponds to the S segment, Clades 1-3 dataset, and linkageInformed scenario.

The datasets correspond to the following clade compositions:

    Clade1and3        Clades 1 and 3
    Clade2            Clade 2
    Clade1-3          Clades 1-3
    Clade1-4          Clades 1-4
    inclUnclassified  Clades 1-4 including unclassified sequences

## ANALYSIS SCENARIOS

### allObserved

Analyses using all available observed sequences.

### TwoPerOutbreak

Analyses using two sequences per outbreak.

### linkageInformed

Analyses incorporating linkage-informed sequences.

## MCMC RANDOM SEEDS

The random number seeds used for each MCMC chain are provided in:

    Table BEAST_Segment_Dataset_Scenario_Seed_XML-filename_for_GitHub.csv
