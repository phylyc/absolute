# Creates abs sif from somatic sif
options(error=recover)
suppressMessages(library(readr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(purrrlyr))
suppressMessages(library(tidyr))
suppressMessages(library(GetoptLong))

abs_dir <- "/cga/scarter/rklein/Projects/SDS/workflows/abs_phy"
project_dir <- "/cga/scarter/rklein/Projects/SDS"
msacsOut_dir <- "/cga/scarter/rklein/Projects/SDS/workflows/somatic_analysis_workflow/data/msacs/pon_1019/output"
#force_call_dir <- "/cga/scarter/ncamarda/sds/merge_snvs_and_indels_workflow/run_04.30.18/pair_mafs"
force_call_dir <- "/cga/scarter/rklein/Projects/SDS/workflows/somatic_analysis_workflow/data/merged_snvs_indels"
sif_fn <- "/cga/scarter/rklein/Projects/SDS/data/sifs/SDS_SIF.tsv"
setwd(project_dir)

sif <- read_tsv(sif_fn) %>%
  select(SID = sample_id, individual_id, gender = SEX_GENOTYPE, sample_type) %>%
  mutate(AllelicCapseg_skew = 0.95) %>%
  mutate(alleliccapseg_tsv = file.path(msacsOut_dir, individual_id, "MSACS", "samples", 
                                       str_c(SID, ".panelphased.abs.v1.4.seg"))) %>%
  mutate(combined_snvs_indels = file.path(force_call_dir, individual_id,
                                   str_c(SID, ".maf"))) %>%
  mutate(gender = ifelse(gender == "XY", "Male", "Female")) %>%
  rename(sample_id = SID) %>%
  filter(sample_type == "Tumor") %>%
  filter(sample_id != "MA01-098-8-3-C1") %>%
  filter(!(individual_id %in% c("SDS1-000192", "SDS1-000003"))) %>%
  select(-sample_type)
## 192 only has a single tumor sample so mutation output is a bit broken
## 196 had a sample that was not specified in the SIF, so data is broken
## 003 is the transplant case which will need to be run in tCR mode

temp <- sif$sample_id

sif <- sif %>% select(-sample_id, -individual_id)

sif <- as_data_frame(sif)
rownames(sif) <- temp

write.table(sif, file.path(abs_dir, "SIF.tsv"), quote = F, sep = "\t", row.names = TRUE, col.names = TRUE)


