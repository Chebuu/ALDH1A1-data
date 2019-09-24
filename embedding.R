source('source_https.R')
source_https('https://raw.github.com/Chebuu/mol2vox/mol2vox.R', 'https://raw.github.com/Chebuu/mol2vox/util.R')
###########
# GLOBALS #
###########
# Max_Response
MAX_DIR <- 'data/Max_Response/'
MAX_FILES <- paste0(MAX_DIR, dir(MAX_DIR))
MAX_CSV <- grep('*.csv', MAX_FILES, value=TRUE)
MAX_SDF <- grep('*.sdf', MAX_FILES, value=TRUE)
MAX_TRAIN_CSV <- grep('train.csv', MAX_CSV, value=TRUE)
MAX_TEST_CSV <- grep('test.csv', MAX_CSV, value=TRUE)
MAX_TRAIN_SDF <- grep('train.sdf', MAX_CSV, value=TRUE)
MAX_TEST_SDF <- grep('test.sdf', MAX_CSV, value=TRUE)

####################
# Simple cartesian #
####################
max.train.csv <- read.csv(MAX_TRAIN_CSV)
max.train.xyz <- lapply(max.train.csv$SMILES, smi2xyz)
max.train.xyz <- lapply(max.train.xyz, xyz2numeric)
dimRanges(max.train.xyz[100:1000])

max.test.csv <- read.csv(MAX_TEST_CSV)
max.test.xyz <- lapply(max.test.csv$SMILES, smi2xyz)

################
# Simple voxel #
################

