library("devtools")
library("GenomicFeatures")
library("wiggleplotr")
load_all("../NMDer/")

#Make TxDb from ensembl
txdb_mm <- makeTxDbFromEnsembl("Mus musculus", server="ensembldb.ensembl.org")
saveDb(txdb_mm, "data-raw/txdb/mus_musculus_txdb.sqlite")

#Load TxDb and exptract exons and cdss
txdb_mm = loadDb("data-raw/txdb/mus_musculus_txdb.sqlite")
exons = exonsBy(txdb_mm, by = 'tx', use.names=TRUE)
cds = cdsBy(txdb_mm, by = 'tx', use.names=TRUE)

#Extract PTBP2 transxripts
selected_tx = c('ENSMUST00000029780','ENSMUST00000197833')
ptbp2_data = list(exons = exons[selected_tx], cdss = cds[selected_tx])
devtools::use_data(ptbp2_data, overwrite = TRUE)


