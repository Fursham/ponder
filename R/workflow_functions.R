#' _ workflow: Load and/or import input files
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW. 
#' This function will import and/or load input data for the workflow and
#' will return a list containing GRanges objects for assembled transcript and 
#' reference annotation, and fasta file for genome sequence
#' 
#' @param file filename or /path/to/file of transcript assemblies. File can be of
#' gtf, gff3 or bed format, and its filename should contain the corresponding extension
#' @param gencode reference annotation of coding transcripts for each gene family. This program
#' contain in-built gencode_basic annotation for GRCh38_hg38 and GRCm38_mm10 assemblies. These
#' assemblies can be loaded by parsing "hg38" or "mm10" into this argument. If users wish to 
#' import other assemblies, please provide filename or /path/to/file of reference CDS assemblies.
#' File can be of gtf, gff3 or bed format, and its filename should contain the corresponding extension
#' @param fasta genome sequence in fasta format. This program contain in-built fasta files for 
#' GRCh38_hg38 and GRCm38_mm10 genome. These sequences can be loaded by parsing "hg38" or "mm10" 
#' into this argument. If users wish to import other genomes, please provide filename or 
#' /path/to/file of fasta files.
#'
#' @return
#' A list containing (1) GRanges object for assembled transcripts, (2) GRanges object for 
#' reference CDS assembly, and (3) DNAstring containing genome sequence
#' @export
#' @import rtracklayer
#' @import Biostrings
#' 

prepareInputs <- function(file, gencode, fasta) {
  
  # import assembled transcripts
    infoLog('Importing assembled transcripts...', logf, quiet)
    fileformat = tail(unlist(strsplit(file, '\\.')), "1")
    if (!fileformat%in%c('gff3','gff','gtf','bed')) {
      stopLog('Incorrect file extension format (Transcript assembly)')
    }
    inputGRanges = rtracklayer::import(file, format = fileformat)

  # load gencode basic GRanges or import from user annotation
  if (any(gencode == c("hg38", "mm10"))) {
    infoLog(sprintf('Loading gencode_basic transcripts from %s assembly...', gencode), logf, quiet)
    
    basicGRanges = gencode_basics[[gencode]]
    fasta = gencode # force pre-loaded fasta to be used if pre-loaded gencode is used
  } else {
    infoLog('Importing user-defined reference annotations...', logf, quiet)
    
    fileformat = tail(unlist(strsplit(gencode, '\\.')), "1")
    if (!fileformat%in%c('gff3','gff','gtf','bed')) {
      stopLog('Incorrect file extension format (Annotation)')
    }
    basicGRanges = rtracklayer::import(gencode, format = fileformat)
  }
  
  # load genome or import user-defined genome sequence
  if (any(fasta == c("hg38", "mm10"))) {
    infoLog(sprintf('Loading genome sequence from %s assembly', fasta), logf, quiet)
    
    genome = genomes[[fasta]]     
  } else {
    infoLog('Importing fasta file...', logf, quiet)
    
    genome = Biostrings::readDNAStringSet(fasta, format="fasta")     
  } 
  return(list(inputGRanges, basicGRanges, genome))
}


#' _ workflow: Testing and correcting transcript assembly
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW. 
#' This function will test and correct the input transcript assembly for consistent chromosome names 
#' and gene ids with CDS assembly and genome sequence.
#' 
#' Non-consistent Gene IDs are corrected by replacing it with Gene IDs from reference CDS assembly
#' that contain overlapping exons. However, this strategy may be inaccurate and may return 
#' false-postive gene ids. To improve the correction, we provide two upper layers of correction
#' that can be initiated by providing the mandatory 'primary_gene_id' and 'secondary_gene_id'
#' arguments
#' 
#' Primary_gene_id refer to the name of the metadata column containing gene ids for each assemnbled transcript. 
#' This is typically 'gene_id'
#' Secondary_gene_id refer to the name of the metadata column that may contain a reference gene id
#' For example: 
#' NEED TO DESCRIBE
#'  
#' But if both args are provided, 
#' function will try to replace 
#' 
#' How to check for metadata columns
#'
#' @param inputGRanges 
#' GRanges object containing transcript and exon features. It is recommended to create this object
#' by importing gtf/gff3 files using rtracklayer::import
#' @param basicGRanges 
#' #' GRanges object containing reference transcript,exon and CDS features. 
#' It is recommended to create this object by importing gtf/gff3 files using rtracklayer::import
#' @param genome
#' DNAstring object containing genome sequence
#' @param primary_gene_id 
#' Optional argument. Typically, this is 'gene_id'. Refer to description for details
#' @param secondary_gene_id 
#' Optional argument. Refer to description for details
#'
#' @return
#' Function will return a modified inputGRanges object
#' 
#' @export
#' @import GenomeInfoDb
#'
#' @examples
preTesting <- function(inputGRanges, basicGRanges, genome, primary_gene_id, secondary_gene_id) {
  
  # testing and correcting chromosome names
  infoLog('Checking and correcting chromosome names...', logf, quiet)
  
  if (seqlevelsStyle(inputGRanges) != seqlevelsStyle(genome)) {
    newStyle <- mapSeqlevels(seqlevels(inputGRanges), seqlevelsStyle(genome))
    newStyle = newStyle[!is.na(newStyle)]
    basicGRanges <- renameSeqlevels(inputGRanges, newStyle)
  }
  if (seqlevelsStyle(basicGRanges) != seqlevelsStyle(genome)) {
    newStyle <- mapSeqlevels(seqlevels(basicGRanges), seqlevelsStyle(genome))
    newStyle = newStyle[!is.na(newStyle)]
    basicGRanges <- renameSeqlevels(basicGRanges, newStyle)
  }

  # testing and correcting gene_ids
  infoLog('Checking and correcting gene_ids...', logf, quiet)
  
  # return error if input and reference do not contain gene_id metadata
  if (is.null(elementMetadata(inputGRanges)$gene_id) | is.null(elementMetadata(basicGRanges)$gene_id)) {
    stop('Gene_ID metadata error: Please ensure input assembled transcripts and/or reference transcript contain gene_id metadata\n')
  }
  
    # preserve gene_ids from input GTF
  elementMetadata(inputGRanges)$old_gene_id <- elementMetadata(inputGRanges)$gene_id
  
  # count number of non standard ID before correction
  nonstand_id_1 = lengths(elementMetadata(inputGRanges) %>% as.data.frame() %>%
                            dplyr::select(gene_id = gene_id) %>%
                            dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
                            dplyr::distinct() %>%
                            as.data.frame())[[1]]
  
  # There are 3 correction steps below.
    # correction 1 and 2 will initiate only if arguments are provided
    # correction 3 will always run, but may give false positive output
    # note that in each step, ONLY the non-standard gene_ids will be corrected
  
  # correction 1: replace primary_gene_id with secondary_gene_id, IF both args are provided
  if (!is.null(primary_gene_id) & !is.null(secondary_gene_id)) {
    infoLog(sprintf('->Attempting to replace %s with %s...', primary_gene_id, secondary_gene_id),
            logf, quiet)
    
    # make dataframe of non_standard primary_gene_ids with its secondary_gene_ids
    gene_id_df = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(primary = primary_gene_id, secondary = secondary_gene_id) %>%
      dplyr::filter(primary%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
      dplyr::distinct() %>%
      as.data.frame()
    
    # replace nonstandard primary_gene_ids from inputGRanges with its secondary_gene_ids
    mcols(inputGRanges)$gene_id <- sapply(mcols(inputGRanges)$gene_id, function(x) {
        ifelse(x%in%gene_id_df$primary, gene_id_df[gene_id_df["primary"] == x][2], x)
      })
  }
  
  # correction 2: replace primary_gene_id with basic gene ID IF:
    # at least primary_gene_id is provided and if it starts with 'ENS'
  if (!is.null(primary_gene_id)) {
    infoLog('-> Attempting to match ensembl gene_ids...', logf, quiet)
    
    # make dataframe of non_standard primary_gene_ids starting with "ENS"
    gene_id_df = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = primary_gene_id) %>%
      dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
      dplyr::distinct() %>%
      dplyr::filter(startsWith(gene_id, "ENS") == TRUE) %>%
      as.data.frame()
    
    # proceed with correction if there is an "ENS" gene_id annotation
    if (lengths(gene_id_df)[[1]] > 0) {
      # make dataframe of annotated gene_id starting with "ENS"
      basic_gene_id_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
        dplyr::select(ens_id = gene_id) %>%
        dplyr::distinct() %>%
        dplyr::filter(startsWith(ens_id, "ENS") == TRUE) %>%
        as.data.frame()
      
      # append suffix of ENS gene_ids on both dataframes 
        # eg. ENSMUST00000029780.11 to ENSMUST00000029780
      gene_id_df$appended_id <-  sapply(strsplit(gene_id_df$gene_id, "\\."), `[`, 1)
      basic_gene_id_df$appended_id <-  sapply(strsplit(basic_gene_id_df$ens_id, "\\."), `[`, 1)

      # combine both data frames based on common 'appended_id' and filter lines with empty 'ens_id'
      gene_id_df = dplyr::left_join(gene_id_df, basic_gene_id_df, by = 'appended_id') %>%
        dplyr::filter(!is.na(ens_id))
      
      # replace gene_id on assembled transcript with ens_id, if match is found
      mcols(inputGRanges)$gene_id = sapply(mcols(inputGRanges)$gene_id, function(x) {
        ifelse(x%in%gene_id_df$gene_id, gene_id_df[gene_id_df["gene_id"] == x][3], x)
      })
    } else {
      warnLog('No ensembl gene ids found\n', logf, quiet)
    }
  }
  
  # correction 3: correct gene_ids by finding overlapping regions.
  infoLog('-> Attempting to correct gene_ids by finding overlapping coordinates...', logf, quiet)
  
  # make dataframe of non_standard gene_ids
  gene_id_df = elementMetadata(inputGRanges) %>% as.data.frame() %>%
    dplyr::select(gene_id = gene_id) %>%
    dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
    dplyr::distinct() %>%
    as.data.frame()
  
  # retrieve first gene_id from gencode that overlaps with first exon of query
  gene_id_df = dplyr::mutate(gene_id_df,
                             ref_gene_id = sapply(gene_id_df$gene_id, function(x) {
                               
                               getfirstOverlap = findOverlaps(
                                 inputGRanges[mcols(inputGRanges)$gene_id == x &
                                                mcols(inputGRanges)$type == "exon"][1],
                                 basicGRanges, 
                                 select = "first")
                               
                               return(ifelse(!is.na(getfirstOverlap), 
                                              mcols(basicGRanges[getfirstOverlap])$gene_id,
                                              NA))
                             })) %>%
    dplyr::filter(!is.na(ref_gene_id))
  
  # replace gene_id on assembled transcript with predicted gene_id from overlap
  mcols(inputGRanges)$gene_id = sapply(mcols(inputGRanges)$gene_id, function(x) {
    ifelse(x%in%gene_id_df$gene, gene_id_df[gene_id_df["gene_id"] == x][2], x)
  }
  )
  
  # count remaining number of non_standard gene_ids and number of corrected ids
  nonstand_id_2 = lengths(elementMetadata(inputGRanges) %>% as.data.frame() %>%
                            dplyr::select(gene_id = gene_id) %>%
                            dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
                            dplyr::distinct() %>%
                            as.data.frame())[[1]]
  
  corrected_ids = nonstand_id_1 - nonstand_id_2
  
  # add a tag to inputGRanges for gene_ids that are standard
  prep_df = elementMetadata(inputGRanges) %>% as.data.frame() %>%
    dplyr::select(gene_id = gene_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(standard_id = ifelse(gene_id%in%unique(mcols(basicGRanges)$gene_id), TRUE, FALSE)) %>%
    as.data.frame()

  mcols(inputGRanges)$standard_id <- sapply(mcols(inputGRanges)$gene_id, function(x) {
    prep_df[prep_df["gene_id"] == x][2]
  })
  
  # report pre-testing analysis and return inputGRanges
  infoLog(sprintf('--> Number of gene_ids corrected: %s', corrected_ids), logf, quiet)
  infoLog(sprintf('--> Remaining number of non-standard gene_ids: %s', nonstand_id_2), logf, quiet)
  if (nonstand_id_2 > 0) {
    warnLog('Transcripts with non-standard gene_ids will skip analysis', logf, quiet)
  }
  return(inputGRanges)
}



#' _ workflow: Preparing for analysis
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will prepare output dataframes and GRanges objects for the _Program 
#'
#' @param inputGRanges 
#' GRanges object containing transcript and exon features. It is recommended to create this object
#' by importing gtf/gff3 files using rtracklayer::import
#' @param basicGRanges 
#' #' GRanges object containing reference transcript,exon and CDS features. 
#' It is recommended to create this object by importing gtf/gff3 files using rtracklayer::import
#'
#' @return
#' Function will return a list containing:
#' (1) dataframe of a list of assembled transcripts and its information for output
#' (2) dataframe of a list of reference CDS transcripts
#' (3) GRangesList of exon coordinates of assembled transcripts grouped by transcript name
#' (4) GRangesList of exon coordinates of reference CDS grouped by transcript name
#' 
#' @export
#' @import dplyr
#' @import GenomicFeatures
#'
#' @examples
prepareAnalysis <- function(inputGRanges, basicGRanges, outdir) {
  
  infoLog('Preparing databases, transcripts and output files...', logf, quiet)
  
  # prepare output dataframe
  report_df = inputGRanges %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, 
                  Original_Gene_ID = old_gene_id,
                  Gene_Name = gene_name,
                  Std_Gene_ID = standard_id,
                  Transcript_ID = transcript_id,
                  Chrom = seqnames,
                  Strand = strand) %>%
    dplyr::mutate(ID_corrected = ifelse(Gene_ID == Original_Gene_ID, FALSE, TRUE),
                  Tx_coordinates = NA,
                  ORF_considered = NA,
                  Alt_tx  = NA,
                  annotatedStart = NA,
                  predictedStart = NA,
                  is_NMD = NA, 
                  dist_to_lastEJ = NA, 
                  uORF = NA,
                  threeUTR = NA) %>%
    dplyr::mutate(CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, AL = NA, IR = NA, 
                   ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, al = NA, ir = NA) %>%
    as.data.frame()
  
  # prepare df of gencode_basic transcripts
  basicTX_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_TX_ID = transcript_id) %>%
    dplyr::distinct() %>% as.data.frame()
  
  combined_report_df = dplyr::right_join(basicTX_df, report_df, by = 'Gene_ID')
  
  # prepare databases
  inputDB = makeTxDbFromGRanges(inputGRanges)
  basicDB = makeTxDbFromGRanges(basicGRanges)
  
  # prepare exon structures by transcripts/CDSs
  inputExonsbyTx = exonsBy(inputDB, by="tx", use.names=TRUE)
  basicExonsbyCDS = cdsBy(basicDB, by="tx", use.names=TRUE)
  basicExonsbyTx = exonsBy(basicDB, by="tx", use.names=TRUE)
  
  
  return(list(combined_report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx))
}






#' _ workflow: Test transcripts for NMD feature
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will analyze the assembled transcripts for NMD features
#'
#' @param report_df 
#' dataframe of a list of assembled transcripts. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param basicTX_df 
#' dataframe of a list of reference CDS transcripts. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param inputExonsbyTx 
#' GRangesList of exon coordinates of assembled transcripts grouped by transcript name. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param basicExonsbyCDS 
#' GRangesList of exon coordinates of reference CDS grouped by transcript name. It is recommended to create this object using
#' prepareAnalysis function from the workflow. see ?prepareAnalysis
#' @param genome 
#' DNAstring containing genome sequence
#' @param PTC_dist 
#' Numerical value referring to minimium distance of premature stop codon to last exon junction
#' @param testNonClassicalNMD 
#' Boolean value (TRUE/FALSE) on whether to test for non-classical NMD features
#'
#' @return
#' updated form of report_df
#' @export
#' @import dplyr
#' @import GenomicRanges
#'
#' @examples
testNMDfeatures <- function(report_df, inputExonsbyTx, basicExonsbyCDS, basicExonsbyTx, genome, PTC_dist = 50, testNonClassicalNMD = FALSE) {
  
  infoLog('Detecting NMD features...', logf, quiet)
  # set a progress bar for analysis
  progressbar = invisible(txtProgressBar(min = 0, max = nrow(report_df), width = 40, style = 3))
  
  
  # test for NMD features for each assembled transcript
  for (i in 1:nrow(report_df)) {
    
    # update progressbar
    if (quiet == FALSE) {
      setTxtProgressBar(progressbar, i); Sys.sleep(0.1); #cat(sprintf('\tTranscripts processed: %s', i))
    }
    
    # add transcript coordinates into report_df as a string of coordinates
    report_df$Tx_coordinates[i] = paste(ranges(inputExonsbyTx[[report_df$Transcript_ID[i]]]),
                                        collapse = ';')
    
    # only analyze transctips with standard gene id
    if (report_df$Std_Gene_ID[i] == TRUE) {
      
      # prepare output list 
      out = list(ORF_considered = NA, 
                 Alt_tx  = NA, 
                 is_NMD = NA, 
                 dist_to_lastEJ = NA, 
                 annotatedStart = NA,
                 predictedStart = NA,
                 uORF= NA, 
                 threeUTR = NA)
      
      # get gencode_basic transcript IDs that has a CDS sequence
        # there could be more than 1 CDS transcripts for each gene
        # transcripts which are non-coding are omitted
      thisbasicCDS = basicExonsbyCDS[[report_df$Ref_TX_ID[i]]]
      
      if (is.null(thisbasicCDS[1])) {
        next
      } else {
        
        infoLog(sprintf('Query : %s  Ref : %s', report_df$Transcript_ID[i], report_df$Ref_TX_ID[i]), logf, quiet = TRUE)
        
        # run test and update output list
        NMDreport = testNMDvsThisCDS(thisbasicCDS, 
                                     inputExonsbyTx[[report_df$Transcript_ID[i]]], 
                                     refsequence = genome, 
                                     PTC_dist, 
                                     testNonClassicalNMD)
        out = modifyList(out, NMDreport)
        
        # update analyzed ORF coordinates into output
        if (!is.na(out$ORF_considered[1])) {
          out$ORF_considered = paste(ranges(out$ORF_considered), collapse = ';')
        }
        
        # classify alternative splicing events
        altevents = getASevents(basicExonsbyTx[[report_df$Ref_TX_ID[i]]], inputExonsbyTx[[report_df$Transcript_ID[i]]])
        if (!is.null(altevents)) {
          report_df[i,] = modifyList(report_df[i,], altevents)
        }
        
        # update transcript information with NMD data
        report_df[i,] = modifyList(report_df[i,], out)
        
      }
    }

  }
  cat('\n')
  return(report_df)
}



#' _ workflow: Test transcript for NMD feature against a reference CDS
#' 
#' @description 
#' THIS FUNCTION IS PART OF THE _ PROGRAM WORKFLOW.
#' This function will compare the query transcript with a reference CDS and tests
#' if insertion of alternative segments into the CDS generates an NMD substrate
#'
#' @param knownCDS 
#' @param queryTx 
#' @param refsequence 
#' @param PTC_dist 
#' @param nonClassicalNMD 
#'
#' @return
#' @export
#'
#' @examples
testNMDvsThisCDS <- function(knownCDS, queryTx, refsequence, PTC_dist = 50, nonClassicalNMD = FALSE) {
  
  # prep output list
  output = list(ORF_considered = NA,
             Alt_tx  = NA,
             is_NMD = NA, 
             dist_to_lastEJ = NA, 
             annotatedStart = NA,
             predictedStart = FALSE,
             uORF = NA,
             threeUTR = NA
             )
  
  # precheck for annotated start codon on query transcript and update output
  pre_report = testTXforStart(knownCDS, queryTx, full.output=TRUE)
  output = modifyList(output, pre_report["annotatedStart"])
  
  # return if there is no shared exons between transcript and CDS
  if (is.na(pre_report$txrevise_out[1])) {
    return(output)
  }
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if ((pre_report$annotatedStart == FALSE) |
      (pre_report$annotatedStart == TRUE & pre_report$firstexlen < 3)) {
    pre_report = reconstructCDSstart(knownCDS = knownCDS, 
                                     queryTx = queryTx,
                                     txrevise_out = pre_report$txrevise_out,
                                     refsequence = refsequence,
                                     full.output = TRUE)

    # return if CDS with new 5' do not contain a start codon
    if (is.na(pre_report$ORF[1])) {
      return(output)
    }
  } 

  # reconstruct CDS with insertion of alternative segments
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out)
  output$Alt_tx = augmentedCDS$Alt_tx
  
  NMDreport = testClassicalNMD(augmentedCDS$ORF, fasta = refsequence, minmPTCdist = PTC_dist, full.output = TRUE)
  
  # update output only if ORF has a stop
  if (length(NMDreport$stopcodons) > 0) {
    output = modifyList(output, list(is_NMD = elementMetadata(NMDreport$stopcodons[1])$is_NMD, 
                                     dist_to_lastEJ = elementMetadata(NMDreport$stopcodons[1])$lastEJ_dist,
                                     ORF_considered = NMDreport$ORF))
    
    # test non Classical NMD only if an ORF have been generated
    if (nonClassicalNMD == TRUE) {
      nonClassicalNMDreport = testNonClassicalNMD(NMDreport$ORF, queryTx, refsequence)
      output = modifyList(output, nonClassicalNMDreport)
    }
  }
  return(output)
}



getASevents <- function(transcript1, transcript2) {
  
  ASlist = list(CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, AL = NA, IR = NA, 
             ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, al = NA, ir = NA)
  
  ASoutput = classifyAltSegments(transcript1, transcript2)
  if (is.null(ASoutput)){
    return(NULL)
  } else if (length(ASoutput[[1]]) == 0 & length(ASoutput[[2]]) == 0) {
    return(NULL)
  }
  
  combinedASoutput = unlist(append(ASoutput[[1]], ASoutput[[2]]))
  for (i in 1:length(combinedASoutput)) {
    eventindex = match(elementMetadata(combinedASoutput[i])$AS_class, names(ASlist))
    if (is.na(ASlist[[eventindex]])) {
      ASlist[[eventindex]] = paste(ranges(combinedASoutput[i]))
    } else {
      ASlist[[eventindex]] = paste(c(ASlist[[eventindex]], paste(ranges(combinedASoutput[i]))), collapse = ';')
    }
  }
  return(ASlist[!is.na(ASlist)])
}





stopLog <- function(text, file) {
  write(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
              sprintf("[STOP] %s", text)), file)
  message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                sprintf("[ERROR] %s", text)))
  stop()
}

infoLog <- function(text, file, quiet = FALSE) {
  
  write(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
              sprintf("[INFO] %s", text)), file)
  if (quiet == FALSE) {
    message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                  sprintf("[INFO] %s", text)))
  }
  Sys.sleep(0.3)
}

warnLog <- function(text, file, quiet = FALSE) {
  
  write(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
              sprintf("[WARN] %s", text)), file)
  if (quiet == FALSE) {
    message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                  sprintf("[WARN] %s", text)))
  }
  Sys.sleep(0.3)
}

