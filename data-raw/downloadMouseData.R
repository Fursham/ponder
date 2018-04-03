library("devtools")
library("GenomicFeatures")
library("wiggleplotr")
library("AnnotationHub")
load_all("../NMDer/")

# load mouse genome sequence
AnnHub <- AnnotationHub()
# query(AnnHub, c("Mus musculus", "release-91"))
mmus_dna <- AnnHub[["AH60492"]] # TwoBitFile for GRCm38 primary assembly
## is there a way to save the twobit file??
devtools::use_data(mmus_dna, overwrite = TRUE)

# Make TxDb from ensembl
txdb_mm <- makeTxDbFromEnsembl("Mus musculus", server="ensembldb.ensembl.org")
saveDb(txdb_mm, "data-raw/txdb/mus_musculus_txdb.sqlite")

# Load TxDb and exptract exons and cdss
txdb_mm = loadDb("data-raw/txdb/mus_musculus_txdb.sqlite")
exons = exonsBy(txdb_mm, by = 'tx', use.names=TRUE)
cds = cdsBy(txdb_mm, by = 'tx', use.names=TRUE)

# Extract Ptbp2 transcripts
ptbp2_testTx = c('ENSMUST00000029780', 'ENSMUST00000197833')
ptbp2_testData = list(exons = exons[ptbp2_testTx], 
                      cdss = cds[ptbp2_testTx], 
                      noNMD = cds$ENSMUST00000029780, 
                      NMD = exons$ENSMUST00000197833,
                      diffstart = exons$ENSMUST00000198399)
devtools::use_data(ptbp2_testData, overwrite = TRUE)

# Extract Bak1 transcripts
bak1_testTx = c('ENSMUST00000078691', 'ENSMUST00000025034') 
bak1_NMDTx = sort(unlist(append(
  reduce(cds$ENSMUST00000078691), 
  reduce(exons$ENSMUST00000025034[5]))), 
  decreasing=TRUE)

bak1_testData = list(exons = exons[bak1_testTx], 
                     cdss = cds[bak1_testTx], 
                     noNMD = cds$ENSMUST00000078691, 
                     NMD = exons$ENSMUST00000025034)
devtools::use_data(bak1_testData, overwrite = TRUE)

# Extract Psd95 transcripts
psd95_testTx = c('ENSMUST00000108589', 'ENSMUST00000123687')
psd95_NMDTx = cds$ENSMUST00000108589[cds$ENSMUST00000108589 != cds$ENSMUST00000108589[20]]
psd95_testData = list(exons = exons[psd95_testTx], 
                      cdss = cds[psd95_testTx], 
                      noNMD = cds$ENSMUST00000108589, 
                      NMD = psd95_NMDTx)
devtools::use_data(psd95_testData, overwrite = TRUE)




