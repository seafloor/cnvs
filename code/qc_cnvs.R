# header ------------------------------------------------------------------
# September 2017
# Matt Smith
# Cardiff University
#
# *outline:*
#   purpose: batch qc for affymetrix arrays
#   
#   input: 
#     - axiom report file (from APT)
#     - .qcsum file (from PennCNV)
#     - annotated rawcnv file (from Elliot's script)
#     - all should be in same folder
#   output: 
#     - single pdf report of batch qc
#     - csv file for importing into spss (needs 'pathogenic_checked' column added)
#
# scripts (try to) follow Google's R style guide: https://google.github.io/styleguide/Rguide.xml


# setup -------------------------------------------------------------------
library(dplyr)
library(skimr)
library(ggplot2)
library(stringr)
library(readr)
library(hexbin)
library(gridExtra)
#library(gridBase)
setwd("~/Documents/phd/biobank_analyses/cnvs")


# import ------------------------------------------------------------------
axiom.report <- read_tsv("AxiomGT1.report.txt", skip = 4780)
qcsum <- read_tsv("Batch1.qcsum")
annotated <- read_tsv("batch1_cnv_annotated.txt")


# reformat ----------------------------------------------------------------
axiom.report <- axiom.report %>%
  mutate(File = str_replace(cel_files, ".CEL", "")) %>%
  select(File, call_rate, computed_gender)

qcsum <- qcsum %>%
  mutate(File = str_replace(File, "split1.", "")) %>%
  mutate(ExcessCNVS = as.factor(ifelse(NumCNV > 10, ">10", "<10"))) %>%
  left_join(axiom.report, by = "File") %>%
  select(File, LRR_SD, WF, NumCNV, call_rate, computed_gender, ExcessCNVS) %>%
  rename(ID = File)

n.people.beforeqc <- nrow(qcsum)
n.cnvs.beforeqc <- sum(qcsum$NumCNV)

# filter people



annotated <- annotated %>%
  mutate(ID = str_replace(ID, "split1.", "")) %>%
  left_join(qcsum, by = "ID") %>%
  mutate(Density = Size / Probe)



# create graphics ---------------------------------------------------------
# visualise output before filtering
by.cnvnum <- qcsum %>%
  select(-File) %>%
  group_by(ExcessCNVS) %>%
  summarise(n = n(),
            mean_WF = mean(WF),
            sd_WF = sd(WF),
            mean_LRRSD = mean(LRR_SD),
            sd_LRRSD = sd(LRR_SD))

wf.heat <- ggplot(data = qcsum, mapping = aes(x = WF, y = NumCNV)) +
  #geom_point(alpha = 0.6) +
  geom_hex(bins = 50) +
  ggtitle("Wave Factor vs. Number of CNVs") +
  xlab("Wave Factor") +
  ylab("Number of CNVs") +
  geom_vline(xintercept = 0.02, col = "red") +
  geom_vline(xintercept = -0.02, col = "red") +
  scale_colour_discrete(name="Number of CNVs") +
  theme_minimal()

lrr.heat <- ggplot(data = qcsum, mapping = aes(x = LRR_SD, y = NumCNV)) +
  #geom_point(alpha = 0.6) +
  geom_hex(bins = 50) +
  ggtitle("s.d. of LRR vs. Number of CNVs") +
  xlab("Standard Deviation of Log R Ratio") +
  ylab("Number of CNVs") +
  geom_vline(xintercept = 0.2, col = "red") +
  scale_colour_discrete(name="Number of CNVs") +
  theme_minimal()

# adding commonly-viewed histograms
hist1 <- ggplot(data = qcsum, mapping = aes(x = WF)) +
  geom_histogram() +
  xlab("Wave Factor") +
  ylab("Count") +
  theme_minimal()

hist2 <- ggplot(data = qcsum, mapping = aes(x = call_rate)) +
  geom_histogram() +
  xlab("Genotyping Call Rate") +
  ylab("Count") +
  theme_minimal()

hist3 <- ggplot(data = qcsum, mapping = aes(x = LRR_SD)) +
  geom_histogram() +
  xlab("Standard Deviation of Log R Ratio") +
  ylab("Count") +
  theme_minimal()

hist4 <- ggplot(data = qcsum, mapping = aes(x = NumCNV)) +
  geom_histogram() +
  xlab("Number of CNVs") +
  ylab("Count") +
  theme_minimal()


# filtering ---------------------------------------------------------------
# filter people
low.callrate <- annotated %>%
  filter(call_rate <= 0.96)

high.numcnvs <- annotated %>%
  filter(NumCNV > 30)

extreme.wf <- annotated %>%
  filter(abs(WF) >= 0.03)

removed.samples <- rbind(low.callrate, high.numcnvs, extreme.wf)
removed.samples <- removed.samples %>%
  filter(ID %in% unique(ID))

# post-people filtering
annotated <- annotated %>%
  filter(call_rate <= 0.96) %>%
  filter(NumCNV > 30) %>%
  filter(abs(WF) >= 0.03)

# filter cnvs

# create pdf --------------------------------------------------------------
# set up pdf layout
pdf.name = "CNV_qc.pdf" # change this name
pdf(file = pdf.name, width = 8.27, height = 11.69) 
par(mfrow = c(6, 1),
    mar = c(2, 5, 2, 1),
    oma = c(3, 0, 4, 0))

# add headings, table and graphics
grid.arrange(wf.heat, lrr.heat, hist1, hist2, hist3, hist4, ncol = 2)
#grid.table(by.cnvnum)

# close pdf
#mtext(paste("\nCNV traces generated using", batch.dir, "on", format(Sys.time(), "%d-%m-%Y"), "using R version", getRversion(), sep = " "), outer = TRUE, side = 1, cex = 0.6)
dev.off() 