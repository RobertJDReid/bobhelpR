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

#' generate a vector of numbers based on start and end points of sequences
#'
#' @param start_vector vector of start numbers for generation of sequences
#' @param end_vector vector of end numbers for generation of sequences
#' @export

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

#' read BLAST results file from the Y1000 BLAST server
#'
#' @param file path to BLAST results text file
#' @param skip number of lines to skip to get to header
#' @param prefix leading prefix text to remove from header line
#' @export

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

#' Find the maximum alignment score on each subject DNA contig contained in the BLAST data
#'
#' @param contig_names vector of contig names taken from BLAST data
#' @param align_mid vector of alignment midpoints from BLAST data
#' @param score vector of BLAST scores from BLAST data
#' @export

maxInContig <- function(contig_names, align_mid, score) {

  # capture calling variable name for contig_names
  #x=last(substitute(contig_names))
  x = as.character(substitute(contig_names))
  y = last(unlist(str_split(x,"\\$")))

  #message(paste0("input contig names are ",x))

  max_positions <- data.frame(contig_names,align_mid, score) %>%
    group_by(contig_names) %>%
    filter(score == max(score)) %>%
    select(contig_names,max_loc = align_mid) %>%
    # restore original column name from calling variable
    #rename(!!x := contig_names)
    rename(!!y := contig_names)
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

#' group fragmented alignments based on proximity to highest score on each subject template
#'
#' @param BLASTdf dataframe of BLAST data from Y1000 project
#' @param query_length_DNA DNA length of the query sequence
#' @export

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

#' Sum scores of grouped alignments
#'
#' @param BLASTdf dataframe containing BLAST results from Y1000 server
#' @export

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

#' find maximum score of alignment groups for each species
#'
#' @param BLASTdf data frame of Y1000 BLAST server information
#' @export

maxGroupScore <- function(BLASTdf) {
  bySpeciesMax <- sumGroupScores(BLASTdf) %>%
    mutate(copy_sub = subject_id) %>%
    separate(copy_sub,c("species",NA),sep=".fas") %>%
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

#' remove sequence insertions from sequence 2 based on gaps in sequence 1
#'
#' @param seq1 string containing a reference sequence
#' @param seq2 string containing a sequence with insertions to remove
#' @export

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

#' Return number segments of aligned sequences taking gaps into account
#'
#' @param start integer start position of the alignment with respect to the query strand
#' @param query string containing query strand sequence from alignment
#' @param subject string containing subject strand sequence from alignment
#' @export

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
# note that this really is broken. If a gap of greater than 1 is made, it ignores anything less that that gap number and does not decrement
# the sequence appropriately
# may be able to fix with addition of negative look ahead and behind in the regex
# for instance this removes only the specified number of repeated 'p's in a string
#  str_remove_all(fruits, "(?<!p)p{3}(?!p)")
# this one removes a range
# str_remove_all(fruits, "(?<!p)p{1,3}(?!p)")

#' find all positions of gaps in a sequence
#'
#' returns a vector of positions in the sequence string representing sequence position after which a gap of \code{gapSize} or greater occur
#' numbers returned are for a completely ungapped sequence
#'
#' @param sequence a sequence string with gaps denoted by '-' from a sequence alignment of BLAST search
#' @param gapSize integer defining the number of gaps >= gapSize to consider when finding positions
#' @importFrom stringr str_locate_all str_remove_all
#' @export

gapPositions <- function(sequence, gapSize = 1) {
  if(gapSize >1) {
    underPattern <- paste0('(?<!-)-{1,',gapSize-1,'}(?!-)')
    #message(underPattern)
    sequence <- str_remove_all(sequence, underPattern)
  }
  pattern <- paste0('(?<!-)-{',gapSize,',}(?!-)')
  positionDF <- as.data.frame(str_locate_all(sequence,pattern)[[1]]) %>%
    mutate(len = end - start + 1,
           cum = cumsum(len),
           recount = end - cum)
  return(positionDF$recount)
}

# gapPositions <- function(sequence, gapSize = 1) {
#   pattern <- paste0('-{',gapSize,',}')
#   positionDF <- as.data.frame(str_locate_all(sequence,pattern)[[1]]) %>%
#     mutate(len = end - start + 1,
#            cum = cumsum(len),
#            recount = end - cum)
#   return(positionDF$recount)
# }

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

#' gets yeast UniProt info from UniProt URL and saves to local file
#'
#' This function simply loads a list of the yeast genes in the Swiss UNIPROT database
#' and saves it as a file
#'
#' @param URL The web location of the file.
#' Defaults to \code{'https://www.uniprot.org/docs/yeast.txt'}
#' @param path The directory path to store the file
#' Defaults to \code{"uniprot.txt"} in working directory.
#' @importFrom RCurl getURL
#' @export

getYeastUniProtInfo <- function (URL='https://www.uniprot.org/docs/yeast.txt',path="uniprot.txt") {
  getURL(URL) %>% readr::write_file(path)
}


#
# ----------------------------------------------------------
#
#                 loads Yeast UniProt info from local file
#
# INPUT
#   path - defaults to "uniprot.txt"
#   skip - defaults to 58
#   colwidths - defaults to c(75, 20, 11, 12, 11, 5, 4, 3)
#   names - defaults to c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")
#
# OUTPUT
#   df containing file info

#' loads Yeast UniProt info from local file
#'
#' Function reads file \code{uniprot.txt} from current directory which contains information
#' on yeast proteins in the Swiss UniProt database. \code{uniprot.txt} is fetched from the
#' Swiss uniprot site with command \code{getYeastUniProtInfo}
#'
#' @param path The directory path to read the file
#' Defaults to \code{"uniprot.txt"} in working directory but can be supplied in function call
#' @param skip - defaults to 58
#' @param max - defaults to 6726
#' @param colwidths - defaults to \code{c(75, 20, 11, 12, 11, 5, 4, 3)}
#' @param names - defaults to \code{c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")}
#' @export


loadYeastUniProtInfo <- function (path="uniprot.txt", skip = 58, n_max = 6726, colwidths = c(75, 20, 11, 12, 11, 5, 4, 3),
                                  colnames = c("gene", "ORF", "swiss-prot-AC","swiss-prot-name","SGD-rec", "size-AA", "3D", "chromosome")) {
  read_fwf("uniprot.txt", skip = skip, n_max = n_max, fwf_widths(colwidths, colnames))

}

#
# ----------------------------------------------------------
#
#                 colonyCount
#
# INPUT
#   img - a single plate image
#   labelText - optional
#
# SIDE EFFECT
#   display image with text
#   display color labels for colonies
#
# OUTPUT
#   number of found objects

#' Identifies colonies on a scanned plate and returns number of colonies
#'
#' @param image A single plate image
#' @param labelText Text to display over image
#' @param maxArea Maximum size cutoff to filter out image artefacts
#' @param minArea Minimum size cutoff to filter non-colony objects
#' @param eccentricity Cutoff for colony objects
#' @export

colonyCount <- function(img,labelText="", maxArea = 800, minArea = 25, eccentricity = 0.6) {
  thr <- img < EBImage::otsu(img)  # automatic threshold level detection
  wat <- EBImage::watershed(EBImage::distmap(thr), tolerance = 1, ext = 1) # defaults


  # Statistics for all detected objects
  obj <- screenmill:::object_features(wat)

  # Review distribution of objects to select filtering parameters
  plot(density(obj$area), xlim = c(0, 1000))
  plot(density(na.omit(obj$eccen)))

  # "Colonies"
  colonies <- obj %>%
    filter(
      area < maxArea,
      area > minArea,
      eccen < eccentricity
    )

  # Unlabel non-colony objects
  result <- as.matrix(wat)
  result[which(!result %in% colonies$obj)] <- 0L

  EBImage::display(img, method = "raster")
  text(x = round(nrow(img),0), y = round(ncol(img),0),
       label=labelText,
       col="firebrick", cex=2)
  EBImage::display(colorLabels(result), method = "raster")
  return(nrow(colonies))
}

#' Computes the angle between three points
#'
#' \code{angle} computes the angle between three points
#' @param A Vector containing the xy-cooydinates of point A
#' @param B Vector containing the xy-cooydinates of point B. This point acts as the vertex of angle ABC
#' @param C Vector containing the xy-cooydinates of point C
#' @return Angle between the three points in radians
#' @export

angle<-function(A, B, C, degrees=FALSE){
  vector1 = B - A
  vector2 = C - B
  # atan2 arguments are (y,x)
  ang <- atan2(vector1[2],vector1[1]) - atan2(vector2[2],vector2[1])
  if (degrees==TRUE){
    ang <- ang * 180 / pi
  }
  return(ang)
}
