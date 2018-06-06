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
#' 
#' 
#' 

prepareInputs <- function(queryfile, ref, fasta, in_format, ref_format) {
  
  # import assembled transcripts
    infoLog('Importing assembled transcripts...', logf, quiet)
  
    if (!file.exists(queryfile)){
      stopLog('Input transcript file do not exist', logf)
    }
  
    # use file extension format if provided by user, else try to extract file format from filename
    if (!is.null(in_format)) {
      fileformat = in_format
    } else {
      fileformat = tail(unlist(strsplit(queryfile, '\\.')), "1")
      if (!fileformat%in%c('gff3','gff','gtf','bed')) {
        stopLog('Incorrect input file extension format', logf)
      } 
    }
    inputGRanges = rtracklayer::import(queryfile, format = fileformat)

  # load provided genome_basic assembly, or import user reference annotation
  if (is(ref, 'GenomicRanges')) {
    infoLog(sprintf('Loading gencode_basic transcripts...'), logf, quiet)
    
    basicGRanges = ref
  } else if (is.character(ref)){
    infoLog('Importing user-provided reference annotations...', logf, quiet)
    
    if (!file.exists(ref)){
      stopLog('Reference annotation file do not exist', logf)
    }
    
    if (!is.null(ref_format)) {
      fileformat = ref_format
    } else {
      fileformat = tail(unlist(strsplit(ref, '\\.')), "1")
      if (!fileformat%in%c('gff3','gff','gtf','bed')) {
        stopLog('Incorrect annotation file extension format', logf)
      }
    }
    basicGRanges = rtracklayer::import(ref, format = fileformat)
  }
  
  
  # import genome sequence
    if (is(fasta, 'BSgenome') | typeof(fasta) == 'S4') {
      infoLog('Loading genome sequence', logf, quiet)
      genome = fasta 
    }
    else if (is.character(fasta)) {
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
preTesting <- function(inputGRanges, basicGRanges, genome, correct_chrom, correct_gene_id, primary_gene_id, secondary_gene_id, clusters) {
  
  # testing and correcting chromosome names on query and annotated transcripts
  infoLog('Checking and correcting chromosome names...', logf, quiet)
  
  if (correct_chrom == TRUE) {
    if (seqlevelsStyle(inputGRanges) != seqlevelsStyle(genome)) {
      newStyle <- mapSeqlevels(seqlevels(inputGRanges), seqlevelsStyle(genome))
      newStyle = newStyle[!is.na(newStyle)]
      inputGRanges <- renameSeqlevels(inputGRanges, newStyle)
      
      if (any(!seqlevels(inputGRanges)%in%seqlevels(genome))) {
        seqlevels(inputGRanges, pruning.mode = 'tidy') <- as.vector(newStyle)
        warnLog('Non-standard chromosome IDs in query were removed', logf, quiet)
      }
    }
    if (seqlevelsStyle(basicGRanges) != seqlevelsStyle(genome)) {
      newStyle <- mapSeqlevels(seqlevels(basicGRanges), seqlevelsStyle(genome))
      newStyle = newStyle[!is.na(newStyle)]
      basicGRanges <- renameSeqlevels(basicGRanges, newStyle)
      
      if (any(!seqlevels(basicGRanges)%in%seqlevels(genome))) {
        seqlevels(basicGRanges, pruning.mode = 'tidy') <- as.vector(newStyle)
        warnLog('Non-standard chromosome IDs in reference were removed', logf, quiet)
      }
    }
  } 
  # if user opts out of this correction service, program will test if there are non-standard IDs and return a warning
  # program will continue
  else {
    if (any(!seqlevels(inputGRanges)%in%seqlevels(genome))) {
      warnLog('Non-standard chromosome IDs in query were found', logf, quiet)
    }
    if (any(!seqlevels(basicGRanges)%in%seqlevels(genome))) {
      warnLog('Non-standard chromosome IDs in reference were found. ', logf, quiet)
    }
  }
  
  # correct gene ids if requested by user
  #   else, calculate the number of unmatched gene_ids and return warning
  if (correct_gene_id == TRUE) {
    inputGRanges = matchGeneIDs(inputGRanges, basicGRanges, 
                                primary_gene_id=primary_gene_id, secondary_gene_id=secondary_gene_id, 
                                makefile = FALSE, clusters = clusters)
  } 
  # if user opts out of this service, program will return warning if there are unmatched gene_ids
  #   program will also update inputGRanges to include match_level and old_gene_id metadata
  else {
    
    # get a list of unmatched gene_ids
    nonstand_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
      dplyr::select(gene_id = gene_id) %>%
      dplyr::filter(gene_id%in%unique(mcols(basicGRanges)$gene_id) == FALSE) %>%
      dplyr::distinct() %>%
      dplyr::mutate(status = TRUE)
    
    if (nrow(nonstand_ids) > 0) {
      
      # return warning
      warnLog(sprintf('%s transcripts have unmatched gene_ids and will not be analyzed', nrow(nonstand_ids)), logf, quiet)
      
      # and add match_level metadata to inputGRanges
      unique_ids = elementMetadata(inputGRanges) %>% as.data.frame() %>%
        dplyr::select(gene_id = gene_id) %>% 
        dplyr::mutate(match_level = 0) %>%
        dplyr::left_join(nonstand_ids, by=c('gene_id'='gene_id')) %>%
        dplyr::mutate_at(vars(match_level), funs(ifelse(!is.na(status), 5, .))) %>%
        dplyr::select(gene_id, match_level)
      mcols(inputGRanges)$match_level = unique_ids$match_level
      mcols(inputGRanges)$old_gene_id = mcols(inputGRanges)$gene_id
      
      rm(unique_ids)
    }
    rm(nonstand_ids)
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
  
  # prepare output dataframe by getting transcript information from inputGRanges
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
                  CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, 
                  ATS = NA, AL = NA, APA = NA, IR = NA,
                  ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, 
                  ats = NA, al = NA, apa = NA, ir = NA, NMDcausing = NA)
  
  # prepare df of gencode_basic transcripts
  basicTX_df = elementMetadata(basicGRanges) %>% as.data.frame() %>%
    dplyr::filter(type == 'transcript') %>%
    dplyr::select(Gene_ID = gene_id, Ref_TX_ID = transcript_id) %>%
    dplyr::distinct() %>% as.data.frame()
  
  # join report_df and basicTX_df and generate a full list of query transcripts to be compared with every reference transcript
  # and add in NMDerIDs
  combined_report_df = dplyr::right_join(basicTX_df, report_df, by = 'Gene_ID') %>%
    dplyr::select(NMDer_ID, Gene_ID, Original_Gene_ID, Match_level, 
                  Gene_Name, Transcript_ID, Ref_TX_ID, Chrom, Strand, 
                  Tx_coordinates, annotatedStart, predictedStart, Alt_tx,
                  ORF_considered, is_NMD, dist_to_lastEJ, uORF, threeUTR, uATG, 
                  uATG_frame, Shared_coverage, CE, MX, A5, A3, AF, AL, ATS, APA, IR, 
                  ce, mx, a5, a3, af, al, ats, apa, ir, NMDcausing) %>% 
    dplyr::mutate(rownum = rownames(combined_report_df)) %>%
    rowwise() %>% dplyr::mutate_at(vars(NMDer_ID), 
                                   funs(paste("NMDer", formatC(as.integer(rownum), width=7, flag="0"), sep=""))) %>%
    dplyr::select(-rownum)
  
  # prepare databases
  inputDB = makeTxDbFromGRanges(inputGRanges)
  basicDB = makeTxDbFromGRanges(basicGRanges)
  
  # prepare exon structures by transcripts/CDSs
  inputExonsbyTx = exonsBy(inputDB, by="tx", use.names=TRUE)
  basicExonsbyCDS = cdsBy(basicDB, by="tx", use.names=TRUE)
  basicExonsbyTx = exonsBy(basicDB, by="tx", use.names=TRUE)
  
  # cleanup unused variables
  rm(list = c('inputDB','basicDB'))
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
testNMDfeatures <- function(report_df, inputExonsbyTx, basicExonsbyCDS, 
                            basicExonsbyTx, genome, PTC_dist = 50, 
                            testNonClassicalNMD = FALSE) {
  
  # this is an internal function for testing NMD features on report_df rowwise
  internalfunc = function(x) {
    
    # convert each line into a list so that elements 
    # can be referenced as thisline$_
    thisline = as.list(x)
    
    # Add in coordinates of query transcripts
    #  this part have to be done this way because some exons are 1 nt in length
    #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
    string = paste(ranges(inputExonsbyTx[[thisline$Transcript_ID]]))
    string = sapply(string, function(y) {
      newsstring = ifelse(!grepl('-', y), paste(c(y, '-', y), collapse = ''), y)
      return(newsstring)
    })
    thisline$Tx_coordinates = paste(string, collapse = ';')
    
    # return if query do not correspond to an annotated gene
    if (thisline$Match_level == 5) {
      return(thisline)
    } else {
    
      # return if reference transcript is not a CDS
      thisbasicCDS = basicExonsbyCDS[[thisline$Ref_TX_ID]]
      if (is.null(thisbasicCDS[1])) {
        return(thisline)
      } else {
        
        # run test and update output list
        NMDreport = testNMDvsThisCDS(thisbasicCDS, 
                                     inputExonsbyTx[[thisline$Transcript_ID]], 
                                     genome, 
                                     PTC_dist, 
                                     testNonClassicalNMD)
        thisline = modifyList(thisline, NMDreport)
        
        # classify alternative splicing events
        altevents = getASevents(basicExonsbyTx[[thisline$Ref_TX_ID]], 
                                inputExonsbyTx[[thisline$Transcript_ID]], 
                                thisline$ORF_considered,
                                thisline$is_NMD)
        if (!is.null(altevents)) {
          thisline = modifyList(thisline, altevents)
        }
        
        # update analyzed ORF coordinates into output
        #  this part have to be done this way because some terminal exons have length of 1
        #  and doing the conventional way will output 'xxx' rather than 'xxx-xxx'
        if (!is.na(thisline$ORF_considered[1])) {
          ORFstringlist = sapply(thisline$ORF_considered, function(x) {
            string = paste(ranges(x))
            ifelse(!grepl('-', string), paste(c(string, '-', string), collapse = ''), string)
          })
          thisline$ORF_considered = paste(ORFstringlist, collapse = ';')
        }
        

        thisline[11:41] = as.character(thisline[11:41])
        # update transcript information with NMD data
        return(thisline)
        #report_df[i,] = modifyList(report_df[i,], thisline)
      }
    }
  }
  
  report_df = report_df %>% 
    rowwise() %>% do(data.frame(internalfunc(.))) %>% 
    dplyr::mutate_all(funs(replace(., . == "NA", NA)))
  
  #report_df = rbind(tempdf, do.call(rbind, report_df))
  
  #cat('\n')
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
  
  # test NMD only if there transcript contain an ORF
  if (!is.na(augmentedCDS$ORF_considered)) {
    NMDreport = testNMD(augmentedCDS$ORF_considered, queryTx, fasta = refsequence, distance_stop_EJ = PTC_dist, other_features = nonClassicalNMD)
    output = modifyList(output, NMDreport)
  }
  return(output)
}






getASevents <- function(transcript1, transcript2, orf, is_NMD) {
  
  # prepare list to be returned
  #ASlist = data.frame(CE = NA, MX = NA, A5 = NA, A3 = NA, AF = NA, ATS = NA, AL = NA, APA = NA, IR = NA,
  #             ce = NA, mx = NA, a5 = NA, a3 = NA, af = NA, ats = NA, al = NA, apa = NA, ir = NA,
  #              NMDcausing = NA)
  
  # get AS classifications. transcript 1 is reference and transcript 2 is query in this case
  ASoutput = classifyAltSegments(transcript1, transcript2)
  
  # return if there is no alternative segments between transcripts
  if (is.null(ASoutput)){
    return(NULL)
  } else if (length(ASoutput[[1]]) == 0 & length(ASoutput[[2]]) == 0) {
    return(NULL)
  }
  
  # combine alternative segments and update ASlist with segment coordinates by matching class annotations with the named ASlist
  combinedASoutput = unlist(append(ASoutput[[1]], ASoutput[[2]]))
  combinedASoutput$AS_class = as.character(combinedASoutput$AS_class)
  
  NMDexon = NA
  if (!is.na(is_NMD)) {
    if (is_NMD == TRUE) {
      
      # check if any of the alt segment overlaps with the last coding exon
      altseg_NMD = combinedASoutput[overlapsAny(combinedASoutput, orf[length(orf)])]
      if (length(altseg_NMD) == 1) {
        
        NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
      } 
      # else, check if any of the alt segment overlaps with orf
      else {
        altseg_NMD = combinedASoutput[overlapsAny(combinedASoutput, range(orf))]
        
        # if only 1 segment overlaps, that should be the NMD-causing exon
        if (length(altseg_NMD) == 1) {
          NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
        } 
        # else, we need to test if each segment is in frame with orf
        else if (length(altseg_NMD) > 1) {
          elementMetadata(altseg_NMD)$size = as.data.frame(width(altseg_NMD))
          
          # filter for in-frame segments
          altseg_NMD[elementMetadata(altseg_NMD)$size %% 3 == 0]
          
          # take the upstream most segment as the NMD-causing segment
          if (length(altseg_NMD) > 0) {
            NMDexon = elementMetadata(altseg_NMD)$AS_class[1]
          } 
        }
      }
    }
    
  } 
  
  #########todo
  elementMetadata(combinedASoutput)$val = as.character(paste(ranges(combinedASoutput)))

  prepout = elementMetadata(combinedASoutput) %>% as.data.frame() %>%
    dplyr::group_by(AS_class) %>% dplyr::summarise(vals = paste(val, collapse=";")) %>% as.data.frame()
  ASlist = dplyr::select(prepout, vals) %>% unlist() %>% setNames(prepout[,1]) %>% as.list()

  out = c(ASlist, NMDcausing = NMDexon, Shared_coverage = ASoutput$Shared_coverage)
  return(out)
}


filterdata <- function(df) {
  
  # simplify dataframe by taking NMDer_ID, Transcript_ID and coverages values and each transcript by decreasing coverage
  newdf = dplyr::select(df, NMDer_ID, Transcript_ID, Shared_coverage) %>% as.data.frame() %>%
    dplyr::arrange(Transcript_ID, desc(Shared_coverage)) %>%
    dplyr::distinct(Transcript_ID, .keep_all = TRUE)
  
  # filter df
  outputdf = dplyr::filter(df, NMDer_ID%in%newdf$NMDer_ID == TRUE)
  return(outputdf)
}


outputAnalysis <- function(report_df, filterbycoverage, other_features, make_gtf, output_dir, input, reference, PTC_dist) {
  
  # filter data based 
  if (filterbycoverage == TRUE) {
    report_df = filterdata(report_df)
  }
  
  # prepare splicing summary
  splicingsummarydf = summSplicing(report_df)
  write.table(splicingsummarydf, file = sprintf("%s/NMDer_splicing_summary.txt", output_dir), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # prepare GTF file, if requested
  if (make_gtf == TRUE) {
    generateGTF(report_df, output_dir)
  }
  
  # prepare output file
  infoLog('Saving analysis report...', logf, quiet)
  if (other_features == FALSE) {
    output_df = report_df %>% select(-uORF, -threeUTR, -uATG, -uATG_frame, -group)
  } else {
    output_df = report_df %>% select(-group)
  }
  write(sprintf('# description; Input: %s; PTC_to_EJ: %snt', 
                tail(unlist(strsplit(input, '/')), '1'), PTC_dist), 
        file = sprintf("%s/NMDer_report.txt", output_dir))
  write.table(output_df, file = sprintf("%s/NMDer_report.txt", output_dir), 
              sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
  unlink(logf)
}


summSplicing <- function(df) {
  
  # extract columns from dataframe corresponding to splicing classes
  splicing_df = df %>% as.data.frame() %>% 
    dplyr::select(NMDer_ID, Gene_ID, Gene_Name, Transcript_ID, Ref_TX_ID,
                  Chrom, Strand, is_NMD, CE:NMDcausing) %>% 
    dplyr::mutate_at(vars(CE:ir), funs(ifelse(is.na(.), 0, lengths((strsplit(as.character(.), ';'))))))

  return(splicing_df)
}


#' Title
#'
#' @param df 
#' @param output_dir 
#'
#' @return
#' @export
#'
#' @examples
generateGTF <- function(df, output_dir) {
  
  infoLog('Writing GTF file...', logf, quiet)

  outfile = sprintf("%s/NMDer.gtf", output_dir)
  write ('#', file = outfile)
  gtfdf = apply(df, 1, function(y) {
    
    tempdf = c()
    
    # get line info and extract exon coords
    exon_coords = unlist(strsplit(as.character(y[['Tx_coordinates']]), ';')) %>% as.data.frame() %>% 
      tidyr::separate('.', into = c('start', 'end'), sep = '-') %>%
      dplyr::mutate(exon_num = row_number())
    
    # prepare transcript info and add to tempdf
    txstart = ifelse(y[['Strand']]== '-', exon_coords[nrow(exon_coords),1], exon_coords[1,1])
    txstop = ifelse(y[['Strand']]== '-',exon_coords[1,2], exon_coords[nrow(exon_coords),2])
    
    
    txmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s; is_nmd %s;', 
                         dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]), y[['is_NMD']])
    tx = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'transcript', start = as.integer(txstart), end = as.integer(txstop),
              score = 1000, strand = as.character(y[['Strand']]), frame = '.', meta = txmetainfo)
    tempdf= rbind(tempdf, tx)


    # prepare exon info and add to tempdf
    
    out = apply(exon_coords, 1, function(x) {
      exmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s; exon_number %s;', 
                           dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]), x[['exon_num']])
      ex = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'exon', start = as.integer(x[['start']]), end = as.integer(x[['end']]),
                score = 1000, strand = as.character(y[['Strand']]), frame = '.', meta = exmetainfo)
      return(ex)
      
      #write(paste(ex, collapse = '\t'), file = outfile, append = TRUE)
    })
    out = do.call(rbind, out)
    tempdf = rbind(tempdf, out)
    
    
    # return cds info if transcript has one
    if (!is.na(y[['ORF_considered']])) {
      cds_coords = unlist(strsplit(as.character(y[['ORF_considered']]), ';')) %>% as.data.frame() %>% 
        tidyr::separate('.', into = c('start', 'end'), sep = '-') %>%
        dplyr::mutate(exon_num = row_number())
      
      # extract information about frames of the CDS
      cdsIRanges = IRanges(start = as.integer(cds_coords[,1]), end = as.integer(cds_coords[,2]))
      frames = c(0, head(cumsum(width(cdsIRanges)%%3)%%3, -1))
      cds_coords$frame = frames
      
      # prepare and add start_codon information
      startcodonstart = ifelse(y[['Strand']]== '-', as.integer(cds_coords[1,2]) - 2, cds_coords[1,1])
      startcodonend = ifelse(y[['Strand']]== '-', cds_coords[1,2], as.integer(cds_coords[1,1])+2)
      startmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                              dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]))
      start = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'start_codon', start = startcodonstart, end = startcodonend,
                   score = 1000, strand = as.character(y[['Strand']]), frame = as.character(0), meta = startmetainfo)
      tempdf= rbind(tempdf, start)
      
      #write(paste(start, collapse = '\t'), file = outfile, append = TRUE)
      
      
      
      # prepare and add stop_codon information
      stopcodonstart = ifelse(y[['Strand']]== '-', cds_coords[nrow(cds_coords),1], as.integer(cds_coords[nrow(cds_coords),2]) - 2)
      stopcodonend = ifelse(y[['Strand']]== '-', as.integer(cds_coords[nrow(cds_coords),1])+2, cds_coords[nrow(cds_coords),2])
      
      stopframe = frames[length(frames)]
      stopmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;',
                             dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]))
      stop = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'stop_codon', start = stopcodonstart, end = stopcodonend,
                  score = 1000, strand = as.character(y[['Strand']]), frame = as.character(stopframe), meta = stopmetainfo)
      tempdf= rbind(tempdf, stop)
      #write(paste(stop, collapse = '\t'), file = outfile, append = TRUE)
      
      
      # prepare CDS cooords information
      out = apply(cds_coords, 1, function(x) {
        cdsmetainfo = sprintf('gene_id = %s; transcript_id %s; gene_name %s;', 
                              dQuote(y[['Gene_ID']]), dQuote(y[['NMDer_ID']]), dQuote(y[['Gene_Name']]))
        cds = list(chrom = as.character(y[['Chrom']]), source = 'NMDer', type = 'CDS', start = as.integer(x[['start']]), end = as.integer(x[['end']]),
                   score = 1000, strand = as.character(y[['Strand']]), frame = as.character(x[['frame']]), meta = cdsmetainfo)
        return(cds)

        #write(paste(cds, collapse = '\t'), file = outfile, append = TRUE)
      })
      out = do.call(rbind, out)
      tempdf = rbind(tempdf, out)
      
    }
    
    return(tempdf)
  })
  gtfdf = do.call(rbind, gtfdf)
  write.table(gtfdf, file = outfile, row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE, append = TRUE)
}







## Below are functions to write and possibly print the warnings and information
#   There are packages out there that can do this professionally, but I want the ability
#   for user to decide whether to print these info into console
#   this is where the quiet argument comes in

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

