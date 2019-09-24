# TODO:
#   - Currently all duplicate CIDs are removed with distinct()
#   - - Fix so one of the duplicates is kept

library(dplyr)
library(ChemmineR)
library(ChemmineOB)
library(assertthat)

source('datautil.R')

aid1030 <- read.table(
  'data/AID_1030_datatable_all.csv',
  header=TRUE,
  sep=',', 
  quote='\"',
  numerals='allow.loss',
  stringsAsFactors=TRUE
)

# Remove description rows and select columns
rows.keep <- 6:nrow(aid1030)
cols.keep <- c('PUBCHEM_SID','PUBCHEM_CID', 'Phenotype', 'Potency', 'Efficacy', 'Max_Response')
aid1030 <- aid1030[rows.keep, cols.keep]

# Clean up fields before subsetting
aid1030 <- aid1030 %>%
  # Remove any row with at least 1 NA value (blank observations remain)
  filter_all(any_vars(!is.na(.))) %>% 
  # Convert factors to character, replace blanks with NA
  mutate_all(~ifelse(. == '', NA_character_, as.character(.))) %>%
  # Convert numeral characters to double
  mutate_at(c('Potency','Efficacy','Max_Response'), ~as.numeric(.)) %>%
  # Convert back to factor
  mutate(Phenotype=as.factor(Phenotype))


################
# Max_Response #
################
# Create the Max_Response dataset
aid1030.max <- aid1030 %>%
  # Remove any row with at least 1 NA value
  filter_all(any_vars(!is.na(.))) %>% 
  # Keep phenotype for sampling and validation
  select(PUBCHEM_CID, Phenotype, Max_Response) %>%
  # !!! TODO:: Not a good solution, all duplicates are removed, should keep 1 duplicate
  distinct(PUBCHEM_CID, .keep_all=TRUE) %>%
  # # Scale by zscore
  # mutate(Max_Response=scale(Max_Response)) %>%
  # # Scale by RMS
  # mutate(Max_Response=scale(Max_Response, center=FALSE)) %>%
  # Shift values to 0 minimum
  mutate(Max_Response=shift0(Max_Response))

# Sample equally from each "Phenotype" category
N <- nrow(aid1030.max)
counts <- catCount(aid1030.max$Phenotype)
weights <- c(rep(1, min(counts)), rep(min(counts)/median(counts), median(counts)), rep(min(counts)/max(counts), max(counts)))

aid1030.max.subset <- aid1030.max %>%
  group_by(Phenotype) %>%
  # sample_frac(1, weight=weights, replace=FALSE)
  # sample_frac(median(counts)/N, replace=FALSE)
  sample_n(min(counts), replace=FALSE)

# Get isomeric SMILES for each CID
max.smiles <- getSMILES.CIDs(aid1030.max.subset$PUBCHEM_CID)
aid1030.max.subset$SMILES <- max.smiles

# Create train/test split
aid1030.max.train <- aid1030.max.subset %>%
  group_by(Phenotype) %>%
  sample_frac(0.95, replace=FALSE)

aid1030.max.test <- aid1030.max.subset %>%
  filter(!(PUBCHEM_CID %in% aid1030.max.train$PUBCHEM_CID))

# Sanity check no duplicate samples between train & test sets
assert_that( !any(aid1030.max.test$PUBCHEM_CID %in% aid1030.max.train$PUBCHEM_CID) )

# Write train & test sets to files
# # Convert train SMILES to SDF and write
aid1030.max.train$PUBCHEM_CID
max.train.CID <- aid1030.max.train$PUBCHEM_CID
max.train.SDFset <- smiles2sdf(aid1030.max.train$SMILES)
cid(max.train.SDFset) <- max.train.CID
write.SDF(max.train.SDFset, 'data/Max_Response/train.sdf')
write.table(aid1030.max.train, 'data/Max_Response/train.csv', sep=',', row.names=FALSE)
# # Convert test SMILES to SDF and write
max.test.CID <- aid1030.max.test$PUBCHEM_CID
max.test.SDFset <- smiles2sdf(aid1030.max.test$SMILES)
cid(max.test.SDFset) <- max.test.CID
write.SDF(max.test.SDFset, 'data/Max_Response/test.sdf')
write.table(aid1030.max.test, 'data/Max_Response/test.csv', sep=',', row.names=FALSE)

###########
# Potency #
###########
# Create the Potency dataset
aid1030.pot <- aid1030 %>%
  # Keep phenotype for sampling and validation
  select(PUBCHEM_CID, Phenotype, Potency) %>%
  # !!! TODO:: Not a good solution, all duplicates are removed, should keep 1 duplicate
  distinct(PUBCHEM_CID, .keep_all=TRUE) %>%
  # Remove any row with at least 1 NA value
  filter_all(any_vars(!is.na(.))) %>% 
  # filter(!is.na(Potency)) %>%
  mutate(Potency=shift0(Potency))

# Down-sample inhibitors to balance classes (N 'Inactive' is 0)
idx.activator <- which(aid1030.pot$Phenotype == 'Activator')
idx.inhibitor <- sample(which(aid1030.pot$Phenotype == 'Inhibitor'), 10000)

aid1030.pot.subset <- rbind(
  aid1030.pot[idx.activator,],
  aid1030.pot[idx.inhibitor,]
)

# Get isomeric SMILES for each CID
aid1030.pot.subset <- aid1030.pot.subset %>%
  mutate(SMILES=getSMILES.CIDs(aid1030.pot.subset$PUBCHEM_CID))

# Create train/test split
aid1030.pot.train <- aid1030.pot.subset %>%
  group_by(Phenotype) %>%
  sample_frac(0.95, replace=FALSE)

aid1030.pot.test <- aid1030.pot.subset %>%
  filter(!(PUBCHEM_CID %in% aid1030.pot.train$PUBCHEM_CID))

# Sanity check no duplicate samples between train & test sets
assert_that( !any(aid1030.pot.test$PUBCHEM_CID %in% aid1030.pot.train$PUBCHEM_CID) )

# Write train & test sets to files
# # Convert train SMILES to SDF and write
pot.train.CID <- aid1030.pot.train$PUBCHEM_CID
pot.train.SDFset <- smiles2sdf(aid1030.pot.train$SMILES)
cid(pot.train.SDFset) <- pot.train.CID
write.SDF(pot.train.SDFset, 'data/Potency/train.sdf')
write.table(aid1030.pot.train, 'data/Potency/train.csv', sep=',', row.names=FALSE)
# # Convert test SMILES to SDF and write
pot.test.CID <- aid1030.pot.test$PUBCHEM_CID
pot.test.SDFset <- smiles2sdf(aid1030.pot.test$SMILES)
cid(pot.test.SDF) <- pot.test.CID
write.SDF(pot.test.SDFset, 'data/Potency/test.sdf')
write.table(aid1030.pot.test, 'data/Potency/test.csv', sep=',', row.names=FALSE)


############
# Efficacy #
############
# Create the Efficacy dataset
aid1030.eff <- aid1030 %>%
  # Keep phenotype for sampling and validation
  select(PUBCHEM_CID, Phenotype, Efficacy) %>%
  # !!! TODO:: Not a good solution, all duplicates are removed, should keep 1 duplicate
  distinct(PUBCHEM_CID, .keep_all=TRUE) %>%
  # Remove any row with at least 1 NA value
  filter_all(any_vars(!is.na(.))) %>% 
  # filter(!is.na(Efficacy)) %>%
  mutate(Efficacy=shift0(Efficacy))

# Down-sample inhibitors to balance classes (N 'Inactive' is 0)
idx.activator <- which(aid1030.eff$Phenotype == 'Activator')
idx.inhibitor <- sample(which(aid1030.eff$Phenotype == 'Inhibitor'), 10000)

aid1030.eff.subset <- rbind(
  aid1030.eff[idx.activator,],
  aid1030.eff[idx.inhibitor,]
)

# Get isomeric SMILES for each CID
aid1030.eff.subset <- aid1030.eff.subset %>%
  mutate(SMILES=getSMILES.CIDs(PUBCHEM_CID))

# Create train/test split
aid1030.eff.train <- aid1030.eff.subset %>%
  group_by(Phenotype) %>%
  sample_frac(0.95, replace=FALSE)

aid1030.eff.test <- aid1030.eff.subset %>%
  filter(!(PUBCHEM_CID %in% aid1030.eff.train$PUBCHEM_CID))

# Sanity check no duplicate samples between train & test sets
assert_that( !any(aid1030.eff.test$PUBCHEM_CID %in% aid1030.eff.train$PUBCHEM_CID) )

# Write train & test sets to files
# # Convert train SMILES to SDF and write
eff.train.CID <- aid1030.eff.train$PUBCHEM_CID
eff.train.SDFset <- smiles2sdf(aid1030.eff.train$SMILES)
cid(eff.train.SDFset) <- eff.train.CID
write.SDF(eff.train.SDFset, 'data/Efficacy/train.sdf')
write.table(aid1030.eff.train, 'data/Efficacy/train.csv', sep=',', row.names=FALSE)
# # Convert test SMILES to SDF and write
eff.test.CID <- aid1030.eff.test$PUBCHEM_CID
eff.test.SDFset <- smiles2sdf(aid1030.eff.test$SMILES)
cid(eff.test.SDFset) <- eff.test.CID
write.SDF(eff.test.SDFset, 'data/Efficacy/test.sdf')
write.table(aid1030.eff.test, 'data/Efficacy/test.csv', sep=',', row.names=FALSE)
