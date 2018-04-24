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

prepareInputs <- function(file, gencode, fasta, in_format, ref_format) {
  
  # import assembled transcripts
    infoLog('Importing assembled transcripts...', logf, quiet)
  
    if (!file.exists(file)){
      stopLog('Input transcript file do not exist', logf)
    }
  
    # use file extension format if provided by user, else try to extract file format from filename
    if (!is.null(in_format)) {
      fileformat = in_format
    } else {
      fileformat = tail(unlist(strsplit(file, '\\.')), "1")
      if (!fileformat%in%c('gff3','gff','gtf','bed')) {
        stopLog('Incorrect input file extension format', logf)
      } 
    }
    inputGRanges = rtracklayer::import(file, format = fileformat)

  # load provided genome_basic assembly, or import user reference annotation
  if (any(gencode == c("hg38", "mm10"))) {
    infoLog(sprintf('Loading gencode_basic transcripts from %s assembly...', gencode), logf, quiet)
    
    basicGRanges = gencode_basics[[gencode]]
  } else {
    infoLog('Importing user-defined reference annotations...', logf, quiet)
    
    if (!file.exists(gencode)){
      stopLog('Refernce annotation file do not exist')
    }
    
    if (!is.null(ref_format)) {
      fileformat = ref_format
    } else {
      fileformat = tail(unlist(strsplit(gencode, '\\.')), "1")
      if (!fileformat%in%c('gff3','gff','gtf','bed')) {
        stopLog('Incorrect annotation file extension format')
      }
    }
    basicGRanges = rtracklayer::import(gencode, format = fileformat)
  }
  
  
  # load genome or import user-defined genome sequence
  if (any(fasta == c("hg38", "mm10"))) {
    infoLog(sprintf('Loading genome sequence from %s assembly', fasta), logf, quiet)
    
    genome = genomes[[fasta]]     
  } else {
    infoLog('Importing fasta file...', logf, quiet)
    
    if (!file.exists(fasta)){
      stopLog('Fasta file do not exist')
    } else {
      genome = Biostrings::readDNAStringSet(fasta, format="fasta")
    }
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
preTesting <- function(inputGRanges, basicGRanges, genome, correct_chrom, primary_gene_id, secondary_gene_id) {
  
  # testing and correcting chromosome names
  infoLog('Checking and correcting chromosome names...', logf, quiet)
  
  if (correct_chrom == TRUE) {
    if (seqlevelsStyle(inputGRanges) != seqlevelsStyle(genome)) {
      newStyle <- mapSeqlevels(seqlevels(inputGRanges), seqlevelsStyle(genome))
      newStyle = newStyle[!is.na(newStyle)]
      inputGRanges <- renameSeqlevels(inputGRanges, newStyle)
      
      if (any(!seqlevels(inputGRanges)%in%seqlevels(genome))) {
        seqlevels(inputGRanges, pruning.mode = 'tidy') <- newStyle
        warnLog('Non-standard chromosome IDs in query were removed', logf, quiet)
      }
    }
    if (seqlevelsStyle(basicGRanges) != seqlevelsStyle(genome)) {
      newStyle <- mapSeqlevels(seqlevels(basicGRanges), seqlevelsStyle(genome))
      newStyle = newStyle[!is.na(newStyle)]
      basicGRanges <- renameSeqlevels(basicGRanges, newStyle)
      
      if (any(!seqlevels(basicGRanges)%in%seqlevels(genome))) {
        seqlevels(basicGRanges, pruning.mode = 'tidy') <- newStyle
        warnLog('Non-standard chromosome IDs in reference were removed', logf, quiet)
      }
    }
  } else {
    if (any(!seqlevels(inputGRanges)%in%seqlevels(genome))) {
      stopLog('Non-standard chromosome IDs in query were found')
    }
    if (any(!seqlevels(basicGRanges)%in%seqlevels(genome))) {
      warnLog('Non-standard chromosome IDs in reference were found. ', logf, quiet)
    }
  }
  inputGRanges = matchGeneIDs(inputGRanges, basicGRanges, 
                               primary_gene_id=primary_gene_id, secondary_gene_id=secondary_gene_id, 
                               workflow = TRUE) 
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
                  Match_level = match_level,
                  Transcript_ID = transcript_id,
                  Chrom = seqnames,
                  Strand = strand) %>%
    dplyr::mutate(NMDer_ID = NA,
                  Tx_coordinates = NA,
                  ORF_considered = NA,
                  Alt_tx  = NA,
                  annotatedStart = NA,
                  predictedStart = NA,
                  is_NMD = NA, 
                  dist_to_lastEJ = NA, 
                  uORF = NA,
                  uATG = NA,
                  uATG_frame = NA,
                  threeUTR = NA) %>%
    dplyr::mutate(Shared_coverage = 0, 
                  CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, ATS = NA, AL = NA, APA = NA, IR = NA,
                  ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, ats = NA, al = NA, apa = NA, ir = NA) %>%
    as.data.frame()
  
  # prepare df of gencode_basic transcripts
  basicTX_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_TX_ID = transcript_id) %>%
    dplyr::distinct() %>% as.data.frame()
  
  combined_report_df = dplyr::right_join(basicTX_df, report_df, by = 'Gene_ID') %>%
    dplyr::select(NMDer_ID, Gene_ID, Original_Gene_ID, Match_level, 
                  Gene_Name, Transcript_ID, Ref_TX_ID, Chrom, Strand, 
                  Tx_coordinates, annotatedStart, predictedStart, Alt_tx,
                  ORF_considered, is_NMD, dist_to_lastEJ, uORF, threeUTR, uATG, uATG_frame, 
                  Shared_coverage, CE, MX, A5, A3, AF, ATS, AL, APA, IR, ce, mx, a5, a3, af, ats, al, apa, ir
                  )
  
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
    
    # parse row information into variable
    thisline = report_df[i,] %>% as.list()
    
    # add NMDer_ID and transcript coordinates into report_df as a string of coordinates
    thisline$NMDer_ID = paste("NMDer", formatC(as.integer(i), width=7, flag="0"), sep="")
    thisline$Tx_coordinates = paste(ranges(inputExonsbyTx[[thisline$Transcript_ID]]),
                                        collapse = ';')
    
    # only analyze transctips with standard gene id
    if (as.integer(thisline$Match_level) < 5) {
      
      # get gencode_basic transcript IDs that has a CDS sequence
        # there could be more than 1 CDS transcripts for each gene
        # transcripts which are non-coding are omitted
      thisbasicCDS = basicExonsbyCDS[[thisline$Ref_TX_ID]]
      
      if (is.null(thisbasicCDS[1])) {
        next
      } else {
        
        infoLog(sprintf('Query : %s  Ref : %s', thisline$Transcript_ID, thisline$Ref_TX_ID), logf, quiet = TRUE)
        
        # run test and update output list
        NMDreport = testNMDvsThisCDS(thisbasicCDS, 
                                     inputExonsbyTx[[thisline$Transcript_ID]], 
                                     refsequence = genome, 
                                     PTC_dist, 
                                     testNonClassicalNMD)
        thisline = modifyList(thisline, NMDreport)
        
        # update analyzed ORF coordinates into output
        if (!is.na(thisline$ORF_considered[1])) {
          ORFstringlist = c()
          for (j in 1:length(thisline$ORF_considered)) {
            string = paste(ranges(thisline$ORF_considered[j]))
            if (!grepl('-', string)) {
              string = paste(c(string, '-', string), collapse = '')
            }
            ORFstringlist = c(ORFstringlist, string)
          }
          thisline$ORF_considered = paste(ORFstringlist, collapse = ';')
        }
        
        # classify alternative splicing events
        altevents = getASevents(basicExonsbyTx[[thisline$Ref_TX_ID]], inputExonsbyTx[[thisline$Transcript_ID]])
        if (!is.null(altevents)) {
          thisline = modifyList(thisline, altevents)
        }
        
        # update transcript information with NMD data
        report_df[i,] = modifyList(report_df[i,], thisline)
        
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
             predictedStart = NA,
             uORF = NA,
             threeUTR = NA,
             uATG = NA,
             uATG_frame = NA
             )
  
  # precheck for annotated start codon on query transcript and update output
  pre_report = testTXforStart(queryTx, knownCDS, full.output=TRUE)
  output = modifyList(output, pre_report["annotatedStart"])
  
  # return if there is no shared exons between transcript and CDS
  if (is.na(pre_report$txrevise_out[1])) {
    return(output)
  }
  
  # attempt to reconstruct CDS for transcripts with unannotated start
  if ((pre_report$annotatedStart == FALSE) |
      (pre_report$annotatedStart == TRUE & pre_report$firstexlen < 3)) {
    pre_report = reconstructCDSstart(queryTx, knownCDS,
                                     refsequence,
                                     pre_report$txrevise_out,
                                     full.output = TRUE)
    output = modifyList(output, list(predictedStart = pre_report$predictedStart))

    # return if CDS with new 5' do not contain a start codon
    if (is.na(pre_report$ORF[1])) {
      return(output)
    }
  } 

  # reconstruct CDS with insertion of alternative segments
  augmentedCDS = reconstructCDS(txrevise_out = pre_report$txrevise_out, fasta = refsequence)
  output = modifyList(output, augmentedCDS)
  
  if (!is.na(augmentedCDS$ORF_considered)) {
    NMDreport = testNMD(augmentedCDS$ORF_considered, queryTx, fasta = refsequence, distance_stop_EJ = PTC_dist, other_features = nonClassicalNMD)
    output = modifyList(output, NMDreport)
  }
  return(output)
}



getASevents <- function(transcript1, transcript2) {
  
  ASlist = list(CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, ATS = NA, AL = NA, APA = NA, IR = NA,
                ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, ats = NA, al = NA, apa = NA, ir = NA)
  
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
  out = c(ASlist, Shared_coverage = ASoutput$Shared_coverage)
  return(out)
}


filterdata <- function(df) {
  
  # simplify dataframe by taking NMDer_ID, Transcript_ID and coverages values and each transcript by decreasing coverage
  newdf = dplyr::select(df, NMDer_ID, Transcript_ID, Shared_coverage) %>% as.data.frame() %>%
    dplyr::arrange(Transcript_ID, desc(Shared_coverage))
  
  # get a list of unique Transcript IDs
  txdf = dplyr::select(df, Transcript_ID) %>% as.data.frame() %>%
    dplyr::distinct()
  
  # update txdf with NMDer_IDs with the highest coverage for each transcript
  txdf$NMDer_ID = sapply(txdf[[1]], function(y) {
    newdf[newdf[2] == y][1]
  })
  
  # filter df
  outputdf = dplyr::filter(df, NMDer_ID%in%txdf$NMDer_ID == TRUE)
  return(outputdf)
}


summSplicing <- function(df) {
  
  # extract columns from dataframe corresponding to splicing classes
  splicing_df = dplyr::select(df, NMDer_ID, Gene_ID, Gene_Name, Transcript_ID, Ref_TX_ID,
                              Chrom, Strand, is_NMD, CE:ir) %>% as.data.frame()
  
  for (i in 9:length(splicing_df)) {
    splicing_df[i] = lapply(splicing_df[i], function(x) {
      ifelse(is.na(x), 0, lengths((strsplit(x, ';'))))
    })
  }
  
  return(splicing_df)
}



generateGTF <- function(df, output_dir) {
  
  infoLog('Writing GTF file...', logf, quiet)
  
  gtfoutput = data.frame(chrom = NULL, source = NULL, type = NULL, start = NULL, end = NULL, score = NULL, 
                   strand = NULL, frame = NULL, meta = NULL)
  
  # for loop
  for (i in 1:nrow(df)) {
    
    tempdf = data.frame(chrom = NULL, source = NULL, type = NULL, start = NULL, end = NULL, score = NULL, 
                        strand = NULL, frame = NULL, meta = NULL)
    
    # get line info and extract exon coords
    line = df[i,]
    exon_coords = unlist(strsplit(line$Tx_coordinates, ';')) %>% as.data.frame() %>% 
      tidyr::separate('.', into = c('start', 'end'), sep = '-')
    
    # prepare transcript info
    txstart = ifelse(line$Strand == '-', exon_coords[nrow(exon_coords),1], exon_coords[1,1])
    txstop = ifelse(line$Strand == '-',exon_coords[1,2], exon_coords[nrow(exon_coords),2])
    txmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s; is_nmd %s;', 
                       dQuote(line$Gene_ID), dQuote(line$NMDer_ID), dQuote(line$Gene_Name), line$is_NMD)
    tx = data.frame(chrom = line$Chrom, source = 'NMDer', type = 'transcript', start = as.integer(txstart), end = as.integer(txstop),
                    score = 1000, strand = line$Strand, frame = '.', meta = txmetainfo)
    tempdf = rbind(tempdf, tx)
    
    # return exon info
    for (j in 1:nrow(exon_coords)) {
      exmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s; exon_number %s;', 
                           dQuote(line$Gene_ID), dQuote(line$NMDer_ID), dQuote(line$Gene_Name), j)
      ex = data.frame(chrom = line$Chrom, source = 'NMDer', type = 'exon', start = as.integer(exon_coords[j,1]), end = as.integer(exon_coords[j,2]),
                      score = 1000, strand = line$Strand, frame = '.', meta = exmetainfo)
      tempdf = rbind(tempdf, ex)
    }
    
    # return cds info
    if (!is.na(line$ORF_considered)) {
      cds_coords = unlist(strsplit(line$ORF_considered, ';')) %>% as.data.frame() %>% 
        tidyr::separate('.', into = c('start', 'end'), sep = '-')
      
      cdsIRanges = IRanges(start = as.integer(cds_coords[,1]), end = as.integer(cds_coords[,2]))
      frames = c(0, head(cumsum(width(cdsIRanges)%%3)%%3, -1))
      
      startcodonstart = ifelse(line$Strand == '-', as.integer(cds_coords[1,2]) - 2, cds_coords[1,1])
      startcodonend = ifelse(line$Strand == '-', cds_coords[1,2], as.integer(cds_coords[1,1])+2)
      
      startmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                            dQuote(line$Gene_ID), dQuote(line$NMDer_ID), dQuote(line$Gene_Name))
      start = data.frame(chrom = line$Chrom, source = 'NMDer', type = 'start_codon', start = startcodonstart, end = startcodonend,
                       score = 1000, strand = line$Strand, frame = 0, meta = startmetainfo)
      tempdf = rbind(tempdf, start)
      
      
      
      
      stopcodonstart = ifelse(line$Strand == '-', cds_coords[nrow(cds_coords),1], as.integer(cds_coords[nrow(cds_coords),2]) - 2)
      stopcodonend = ifelse(line$Strand == '-', as.integer(cds_coords[nrow(cds_coords),1])+2, cds_coords[nrow(cds_coords),2])
      stopframe = frames[length(frames)]
      
      stopmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                              dQuote(line$Gene_ID), dQuote(line$NMDer_ID), dQuote(line$Gene_Name))
      stop = data.frame(chrom = line$Chrom, source = 'NMDer', type = 'stop_codon', start = stopcodonstart, end = stopcodonend,
                         score = 1000, strand = line$Strand, frame = stopframe, meta = stopmetainfo)
      tempdf = rbind(tempdf, stop)
      
      
      for (k in 1:nrow(cds_coords)) {
        cdsmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                             dQuote(line$Gene_ID), dQuote(line$NMDer_ID), dQuote(line$Gene_Name))
        cds = data.frame(chrom = line$Chrom, source = 'NMDer', type = 'CDS', start = as.integer(cds_coords[k,1]), end = as.integer(cds_coords[k,2]),
                        score = 1000, strand = line$Strand, frame = frames[k], meta = cdsmetainfo)
        tempdf = rbind(tempdf, cds)
      }
    }
   
    # sort table
    if (line$Strand == '-') {
      tempdf = dplyr::arrange(tempdf, desc(end), type)
    } else {
      tempdf = dplyr::arrange(tempdf, start, type)
    }
    gtfoutput = rbind(gtfoutput, tempdf)
  }
  
  #write table
  write.table(gtfoutput, file = sprintf("%s/NMDer.gtf", output_dir), 
              sep = "\t", col.names = FALSE, row.names = FALSE, 
              qmethod = "double", quote = FALSE)
}









stopLog <- function(text, file) {
  write(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
              sprintf("[STOP] %s", text)), file)
  message(paste(sprintf("[%s]", format(Sys.time(), "%Y-%b-%d %X")),
                sprintf("[ERROR] %s", text)))
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
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

