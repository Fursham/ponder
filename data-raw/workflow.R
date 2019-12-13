library(pondeR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(wiggleplotr)
library(GenomicFeatures)
library(dplyr)

checkInputs(query_gtf, ref_gtf, Mmusculus)

query_gtf <- matchSeqLevels(query_gtf, ref_gtf)
query_gtf<- matchGeneIDs(query_gtf, ref_gtf)

#make exonsby
query_exons = exonsBy(query_gtf, by="tx", use.names=TRUE)
ref_cds = cdsBy(ref_gtf, by="tx", use.names=TRUE)
ref_exons = exonsBy(ref_gtf, by="tx", use.names=TRUE)

#get coverage
query_ids = query_exons %>% as.data.frame() %>%
  select(gene_id, transcript_id) %>%
  distinct()
ref_ids = ref_exons %>% as.data.frame() %>%
  select(gene_id, ref_transcript_id = transcript_id) %>%
  distinct()
q2r <- left_join(query_ids, ref_ids)

q2rcovs <- getCoverages(query_exons, ref_exons, q2r)

#getCDS
query_cds <- predictCDS(query_exons, ref_exons, 
                        Mmusculus, q2rcovs, coverage = 3)

#refine uORF and uATG

#testNMD
predictNMD(query_exons, query_cds)