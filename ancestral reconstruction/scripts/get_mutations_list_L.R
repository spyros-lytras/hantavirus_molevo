#Get reconstructed ancestor sequences for the M segment from IQ-TREE output
#Ancestral nodes inputted from tree manually

library(Biostrings)

setwd("../iqtree_output/")

state1 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state1<- state1[state1$Node=="Node13",]
anc_seq1 <- paste(state1$State, collapse = "")
names(anc_seq1)<- "outbreak_Hondius_anc"

state2 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state2<- state2[state2$Node=="Node12",]
anc_seq2 <- paste(state2$State, collapse = "")
names(anc_seq2)<- "outbreak_Hondius_ref"

state3 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state3<- state3[state3$Node=="Node16",]
anc_seq3 <- paste(state3$State, collapse = "")
names(anc_seq3)<- "outbreak_PV_anc"

state4 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state4<- state4[state4$Node=="Node11",]
anc_seq4 <- paste(state4$State, collapse = "")
names(anc_seq4)<- "outbreak_PV_ref"

state5 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state5<- state5[state5$Node=="Node53",]
anc_seq5 <- paste(state5$State, collapse = "")
names(anc_seq5)<- "outbreak_MN_anc"

state6 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state6<- state6[state6$Node=="Node23",]
anc_seq6 <- paste(state6$State, collapse = "")
names(anc_seq6)<- "outbreak_MN_ref"

state7 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state7<- state7[state7$Node=="Node27",]
anc_seq7 <- paste(state7$State, collapse = "")
names(anc_seq7)<- "outbreak_epyugen_anc"

state8 <- read.table("./sequences_L_ex_lab_cells_2026_05_27v7.cds.linsi.fasta.state", header = TRUE, stringsAsFactors = FALSE)
state8<- state8[state8$Node=="Node26",]
anc_seq8 <- paste(state8$State, collapse = "")
names(anc_seq8)<- "outbreak_epyugen_ref"

seqs <- c(
  anc_seq1,
  anc_seq2,
  anc_seq3,
  anc_seq4,
  anc_seq5,
  anc_seq6,
  anc_seq7,
  anc_seq8
)

writeLines(
  unlist(
    lapply(names(seqs), function(nm) {
      c(paste0(">", nm), seqs[[nm]])
    })
  ),
  "../ancestral_sequences/ancestral_sequences_L.fasta"
)

# Read fasta containing ancestral sequences aligned to reference sequence of Orthohantavirus Andes NC_003468.2
fasta <- readDNAStringSet("../ancestral_sequences_aligned/ancestral_sequences_L_aligned_to_RefSeq.fasta")

names(fasta)

anc_seq1 <- as.character(fasta["outbreak_Hondius_anc"])
anc_seq2 <- as.character(fasta["outbreak_Hondius_ref"])
anc_seq3 <- as.character(fasta["outbreak_PV_anc"])
anc_seq4 <- as.character(fasta["outbreak_PV_ref"])
anc_seq5 <- as.character(fasta["outbreak_MN_anc"])
anc_seq6 <- as.character(fasta["outbreak_MN_ref"])
anc_seq7 <- as.character(fasta["outbreak_epyugen_anc"])
anc_seq8 <- as.character(fasta["outbreak_epyugen_ref"])


#Codon start and end positions assigned with respect to reference
syn_nonsyn_table <- function(ref_seq, query_seq, cds_start = 36, cds_end = 6497) {
  ref_seq <- toupper(ref_seq)
  query_seq <- toupper(query_seq)
  
  if (nchar(ref_seq) != nchar(query_seq)) {
    stop("Sequences must have the same length.")
  }
  if (cds_start < 1 || cds_end > nchar(ref_seq) || cds_start > cds_end) {
    stop("Invalid cds_start/cds_end coordinates.")
  }
  
  # Extract CDS region
  ref_cds <- substring(ref_seq, cds_start, cds_end)
  qry_cds <- substring(query_seq, cds_start, cds_end)
  
  # Trim to full codons
  trim_len <- nchar(ref_cds) - (nchar(ref_cds) %% 3)
  ref_cds <- substring(ref_cds, 1, trim_len)
  qry_cds <- substring(qry_cds, 1, trim_len)
  
  codon_starts <- seq(1, trim_len, by = 3)
  ref_codons <- substring(ref_cds, codon_starts, codon_starts + 2)
  qry_codons <- substring(qry_cds, codon_starts, codon_starts + 2)
  
  # Only translate clean codons
  valid <- grepl("^[ACGT]{3}$", ref_codons) & grepl("^[ACGT]{3}$", qry_codons)
  
  ref_aa <- rep(NA_character_, length(ref_codons))
  qry_aa <- rep(NA_character_, length(qry_codons))
  
  if (any(valid)) {
    ref_aa[valid] <- as.character(translate(DNAStringSet(ref_codons[valid]) , genetic.code = GENETIC_CODE,
                                            no.init.codon = TRUE))
    qry_aa[valid] <- as.character(translate(DNAStringSet(qry_codons[valid]) , genetic.code = GENETIC_CODE,
                                            no.init.codon = TRUE))
  }
  
  # Keep only codons that differ
  changed <- ref_codons != qry_codons
  idx <- which(changed)
  
  if (length(idx) == 0) {
    return(data.frame())
  }
  
  # Return one row per differing nucleotide position within a codon
  out_list <- lapply(idx, function(i) {
    ref_bases <- strsplit(ref_codons[i], "")[[1]]
    qry_bases <- strsplit(qry_codons[i], "")[[1]]
    diff_pos <- which(ref_bases != qry_bases)
    
    if (length(diff_pos) == 0) {
      return(NULL)
    }
    
    lapply(diff_pos, function(p) {
      data.frame(
        cds_nt_pos = codon_starts[i] + p - 1,
        genome_nt_pos = cds_start + codon_starts[i] + p - 2,
        codon_num = ((codon_starts[i] - 1) %/% 3) + 1,
        codon_pos = p,
        ref_codon = ref_codons[i],
        qry_codon = qry_codons[i],
        ref_aa = ref_aa[i],
        qry_aa = qry_aa[i],
        effect = if (!valid[i]) {
          "ambiguous"
        } else if (ref_aa[i] == qry_aa[i]) {
          "synonymous"
        } else {
          "nonsynonymous"
        },
        stringsAsFactors = FALSE
      )
    })
  })
  
  out <- do.call(rbind, unlist(out_list, recursive = FALSE))
  
  rownames(out) <- NULL
  out
}

# Get number of mutation per outbreak 
mut1 <- syn_nonsyn_table(anc_seq2, anc_seq1, cds_start = 36, cds_end = 6497)
mut2 <- syn_nonsyn_table(anc_seq4, anc_seq3, cds_start = 36, cds_end = 6497)
mut3 <- syn_nonsyn_table(anc_seq6, anc_seq5, cds_start = 36, cds_end = 6497)
mut4 <- syn_nonsyn_table(anc_seq8, anc_seq7, cds_start = 36, cds_end = 6497)


# Assign outbreak names to each list of mutations
mut1$outbreak <- "Hondius_outbreak"
mut2$outbreak <- "PV_outbreak"
mut3$outbreak <- "MN_outbreak"
mut4$outbreak <- "Epyugen_outbreak"
all_mut <- rbind(mut1, mut2, mut3, mut4)

write.csv(all_mut,"./../mutations list/all_mutations_L_segment.csv", row.names=FALSE)

