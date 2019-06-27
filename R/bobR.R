# First R Package bobhelpR
#
# Just a bunch of helper functions
# but in particular functions to deal with Y1000 BLAST data
#
# ----------------------------------------------------------
#
#                   seqPile function
#
# makes a pileup vector from a vector of sequence starts and a vector of sequence ends
#  as in DNA or protein sequences
#
# Input
#   a vector of start numbers
#   a vector of end numbers
#
# Returns a vector of the collapsed sequence runs for compiling a histogram
#

# returning as a vector of numbers

seqPile <- function(start_vector,end_vector) {
  listPile <- purrr::map2(start_vector,end_vector,seq) %>%
    purrr::reduce(c)
}

#
# ----------------------------------------------------------
#
#                readYKBLAST function
#
# reads a full BLAST report from the Y1000 BLAST server and parses into a dataframe
#
# Input
#   dataDir defaults to 'data/'
#   filename defaults to 'sequenceserver-full_tsv_report.tsv'

readYKBLAST <- function(file, skip = 3, prefix = "# Fields: ") {

    header <- read_lines(file, #pluck header line out of BLAST results file
                             skip = skip, n_max = 1) %>%
    str_remove(prefix) %>% # and clean up the header names
    str_split(", ") %>%
    unlist() %>%
    str_remove_all("\\.") %>%
    str_replace_all(" ","_")

  # BLAST data file

  BLAST_data <- read_delim(file,
                           "\t", escape_double = FALSE, col_names = FALSE,
                           comment = "#", trim_ws = TRUE, na = c("", "NA", "N/A")) %>%
    set_names(header) %>%
    # select(-subject_tax_ids:-subject_super_kingdoms,   # remove columns with no data
    #        -subject_ids, -subject_gis,-subject_accs) %>% # remove seemingly duplicate columns
    separate(`query/sbjct_frames`, c("query_frame","subject_frame"), sep = '/') %>%
    mutate(s_mid = s_start + round((s_end - s_start)/2,digits = 0))
  return(BLAST_data)
}


#
# ----------------------------------------------------------
#
#                      alignmentMidpoints
#
# returns the midpoint position of the alignment wrt the subject template
#

# maxPerContigMidpoint function
#
# calculates the position on a template of the maximum alignment score from a group of contigs and scores
# useful for results from genome-wide alignment searches where there may be multiple "hits" on a contig
#
# INPUT
#   for vectors for variables:
#      contig_names
#      align_start
#      align_end
#      score
#
# OUTPUT
#   df of contig names and the midpoint wrt contig of the max BLAST hit

maxInContig <- function(contig_names, align_mid, score) {

  # capture calling variable name for contig_names
  x=last(substitute(contig_names))

  #message(paste0("input contig names are ",x))

  max_positions <- data.frame(contig_names,align_mid, score) %>%
    group_by(contig_names) %>%
    filter(score == max(score)) %>%
    select(contig_names,max_loc = align_mid) %>%
    # restore original column name from calling variable
    rename(!!x := contig_names)
}

#
# ----------------------------------------------------------
#
#               alignmentGroup function
#
# INPUT
#   A Y1000 BLAST alignment file
# OUTPUT
#   A BLAST alignment file filtered to include max scoring hit on each DNA subject
#   plus lower scoring alignments that are within query DNA length of the high scoring alignment

alignmentGroups <- function (BLASTdf, query_length_DNA) {

  maxPositions <- maxInContig (BLASTdf$subject_id,BLASTdf$s_mid,BLASTdf$score)

  BLASTdf %>% left_join(maxPositions) %>%
    mutate (dif_to_max = abs(s_mid - max_loc)) %>%
    filter(dif_to_max < query_length_DNA)
}

#
# ----------------------------------------------------------
#
#                  scoreSUM function
#
# INPUT
#   A list of DNA subjects - contig, scaffold, chromosome, etc
# OUTPUT
#   DF with summation of grouped BLAST alignment scores

sumGroupScores <- function(BLASTdf) {
  score_sums <- BLASTdf %>%
    group_by(subject_id) %>%
    summarise(
      score_sum = sum(score),
      n = n()
    )
}

#
# ----------------------------------------------------------
#
#                   maxGroupScores

maxGroupScore <- function(BLASTdf) {
  bySpeciesMax <- sumGroupScores(BLASTdf) %>%
    mutate(copy_sub = subject_id) %>%
    separate(copy_sub,c("species",NA),sep=".fas") %>% #somehow I could not get str_split to do this, so I just bashed a column
    group_by(species) %>%
    filter(score_sum == max(score_sum))
}


#
# ----------------------------------------------------------
#
#                       unGap
#
# INPUT
#   seq1
#   seq2
#
# OUTPUT
#   seq2 with insertions wrt seq1 removed

unGap <- function(seq1, seq2) {

  #validate that seq1 has gaps, otherwise return seq2
  if (str_detect(seq1, "-", negate = T)) {
    return (seq2)
  } else {
    char_seq1 <- str_split(seq1,"")[[1]]
    char_seq2 <- str_split(seq2,"")[[1]]
    return (paste0(char_seq2[char_seq1 != "-"], collapse=""))
  }

}

#
# ----------------------------------------------------------
#
#                     fineSegment
#
# INPUT
#   query sequence match start
#   query sequence from match
#   subject sequence from match
#
# Output
#   A vector of positions conserved wrt the query sequence
#

fineSegment <- function(start,query,subject) {
  new.string <- unGap(query,subject)
    # gaps in query now removed from subject
  segments <- as.data.frame(str_locate_all(new.string,'[A-Z]{1,}'))
  pad.segments <- data.frame(
    start = segments$start+start -1,
    end = segments$end + start -1
  )
  segmentList <- seqPile(pad.segments$start, pad.segments$end)

  return(segmentList)
}

#
# ----------------------------------------------------------
#
#                      gapPositions
#
# INPUT
#   sequence with dash character as gap
#
# OUTPUT
#   positions in seqeunce that are opened for gap from - gap opens after numbered position
#

# gapPositions <- function(sequence) {
#   positionDF <- as.data.frame(str_locate_all(sequence,'-{1,}')[[1]]) %>%
#     mutate(len = end - start + 1,
#            cum = cumsum(len),
#            recount = end - cum)
#   return(positionDF$recount)
#
# }

gapPositions <- function(sequence, gapSize = 1) {
  pattern <- paste0('-{',gapSize,',}')
  positionDF <- as.data.frame(str_locate_all(sequence,pattern)[[1]]) %>%
    mutate(len = end - start + 1,
           cum = cumsum(len),
           recount = end - cum)
  return(positionDF$recount)
}

#
# ----------------------------------------------------------
#
#                  getYeastUniProtInfo
#
# INPUT
#   URL - defaults to 'https://www.uniprot.org/docs/yeast.txt'
#   path - defaults to "uniprot.txt"
#   skip - defaults to 58
#   widths - defaults to c(75, 20, 11, 12, 11, 5, 4, 3)
#   names - defaults to c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")
#
# OUTPUT
#   side effect of writing the uniprot.txt file into the working directory

#' get yeast UniProt info
#'
#' This function simply loads a list of the yeast genes in the Swiss UNIPROT database
#' and saves it as a file
#'
#' @param URL The web location of the file.
#' Defaults to \code{'https://www.uniprot.org/docs/yeast.txt'}
#' @param path The directory path to store the file
#' Defaults to \code{"uniprot.txt"} in working directory.
#' @importFrom readr write_file
#' @importFrom RCurl getURL
#' @export

getYeastUniProtInfo <- function (URL='https://www.uniprot.org/docs/yeast.txt',path="uniprot.txt") {
  RCurl::getURL(URL) %>% readr::write_file(path)
}


#
# ----------------------------------------------------------
#
#                 loadYeastUniProtInfo
#
# INPUT
#   path - defaults to "uniprot.txt"
#   skip - defaults to 58
#   colwidths - defaults to c(75, 20, 11, 12, 11, 5, 4, 3)
#   names - defaults to c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")
#
# OUTPUT
#   df containing file info

#' loadYeastUniProtInfo
#'
#' Function reads file \code{uniprot.txt} from current directory which contains information
#' on yeast proteins in the Swiss UniProt database. \code{uniprot.txt} is fetched from the
#' Swiss uniprot site with command \code{getYeastUniProtInfo}
#'
#' #' @param path The directory path to read the file
#' Defaults to \code{"uniprot.txt"} in working directory but can be supplied in function call
#' @param skip - defaults to 58
#' @param colwidths - defaults to \code{c(75, 20, 11, 12, 11, 5, 4, 3)}
#' @param names - defaults to \code{c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")
#' @export


loadYeastUniProtInfo <- function (path="uniprot.txt", skip = 58, colwidths = c(75, 20, 11, 12, 11, 5, 4, 3),
                                  colnames = c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")) {
  read_fwf("uniprot.txt", skip = 58,fwf_widths(colwidths, colnames))

}
