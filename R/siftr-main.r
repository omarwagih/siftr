
if(F){
  require(seqinr)
  require(data.table)
  require(parallel)
  setwd('~/Development/siftr/R')
  source('siftr-helper.R')
}

#' Predict sift scores from a protein alignment
#' 
#' @param aln path to alignment file in fasta format or character vector of protein sequences
#' @param cores number of cores for parallel processing (default: 1)
#' 
#' @return matrix of sift scores
#' 
#' @import data.table
#' @importFrom seqinr read.fasta
#' @importFrom parallel mclapply
#' @export
#' 
#' @examples
#' # Generate dummy alignment 
#' aln = c('SSSS', 'STTT', 'SYSY', 'ASST', 'KKHS')
#' 
#' # Generate matrix of sift scores
#' sift_mat = predictFromAlignment(aln)
#' 
#' # Show summary of run
#' summary(sift_mat)
#' 
#' # Convert to data frame and filter out unreliable predictions
#' sift_df = filterPredictions(sift_mat, score_thresh = 0.05, residue_thresh = 2)
predictFromAlignment <- function(aln, cores=1){
  # Amino acid alphabet
  AA = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
  
  base = system.file("extdata", "precomputed", package = "siftr")
  
  # List of valid amino acids, sorted
  frq = readRDS(file.path(base, 'aa-freq.rds'))
  # Get rank matrix
  rank_matrix = readRDS(file.path(base, 'rank-matrix.rds'))
  # Load precomputed dirichelet data
  diri = readRDS(file.path(base, 'diri.rds'))
  
  
  # Read sequence alignment
  if(is.character(aln) & length(aln) > 1){
    # Alignment is given as character vector
    aln_fa = strsplit(aln, '')
    # We don't have an alignment file
    aln = NA
  }else{
    # Normalize the path
    aln = normalizePath(aln)
    # Assume aln is the path to fasta file
    aln_fa = read.fasta(aln, forceDNAtolower=F, seqtype='AA')
  }
  
  # Number of positions
  npos = checkUnequalSequenceLengths(aln_fa)
  
  # Remove identical seqs
  aln_fa = removeIdenticalSeqs(aln_fa)
  
  
  # Make letter matrix
  aln_mat = letterMatrix(aln_fa, AA, npos)
  
  # Compute median information content for each position (bottleneck)
  mic = medianInfoContent(aln_mat, frq, cores)
  
  # Get pb weights
  pb_weights = pbWeights(aln_mat, return_all = TRUE)
  
  # Use to construct weight matrix
  pwm = makePWM(aln_mat, pb_weights$pb_mat, frq)
  
  # Normalize pb weights
  pb_weights_norm = pbWeightsNorm(pb_weights, aln_mat)
  
  # Get epsilon values, used to construct final matrix
  eps = getEpsilonValues(pb_weights, pb_weights_norm, rank_matrix)
  
  # Make final sift matrix!
  sift_mat = makeSiftMatrix(pb_weights, pb_weights_norm, eps, diri)
  
  rownames(sift_mat) = AA
  colnames(sift_mat) = 1:ncol(sift_mat)
  
  
  attr(sift_mat, 'aln_file') = aln
  attr(sift_mat, 'info_content') = mic
  attr(sift_mat, 'diff_aas') = pb_weights$diff_aas
  attr(sift_mat, 'nseq') = pb_weights$nseq
  attr(sift_mat, 'npos') = pb_weights$npos
  attr(sift_mat, 'num_aas') = colSums( pb_weights$pfm )
  attr(sift_mat, 'query') = aln_fa[[1]]
  attr(sift_mat, 'aa') = AA
  
  class(sift_mat) = c('siftr_matrix', 'matrix')
  sift_mat
}

#' Convert siftr matrix to data frame and filter out unreliable predictions. 
#' For more information on the filtering see \url{http://sift.bii.a-star.edu.sg/www/SIFT_help.html}
#' 
#' @param sift_mat siftr matrix generated using predictFromAlignment
#' @param score_thresh sift score threshold, scores above this threshold are removed (default: 0.05)
#' @param ic_thresh information content threshold, predictions in positions that have an ic above this threshold are considered low confidence and removed (default: 3.25)
#' @param residue_thresh minimum number of aligned residues at the position of prediction (default: 2)
#' 
#' @export
#' @return data frame of predictions
filterPredictions <- function(sift_mat, score_thresh=0.05, ic_thresh=3.25, residue_thresh=2){
  # Amino acid alphabet
  AA = attr(sift_mat, 'aa')
  
  # Number of amino acids / positions
  naa = nrow( sift_mat )
  npos = ncol( sift_mat )
  
  n_aas = attr(sift_mat, 'num_aas')
  n_seq = attr(sift_mat, 'nseq')
  query = attr(sift_mat, 'query')
  median_ic = attr(sift_mat, 'info_content')
  
  dt = data.table(pos=rep(1:npos, each=naa), 
                  ref=rep(query, each=naa), 
                  alt=AA[rep(1:naa, times=npos )], 
                  score=as.double(sift_mat), 
                  median_ic = rep(median_ic, each=naa),
                  n_aa=rep(n_aas, each=naa),
                  n_seq=n_seq)
  dt = dt[dt$score <= score_thresh & dt$n_aa >= residue_thresh & dt$median_ic <= ic_thresh, ]
  as.data.frame(dt)
}

#' Print siftr matrix
#' 
#'  @param x object of class siftr_matrix 
#'  @param list_attributes if TRUE will show stored attributes for object
print.siftr_matrix  <- function (x, list_attributes = FALSE, ...) {
  if(!list_attributes){
    attributes_all <- names(attributes(x))
    attributes_rm <- attributes_all[!(attributes_all %in% c("dim",
                                                            "names", "dimnames", 
                                                            "row.names", "col.names"))]
    attributes(x)[attributes_rm] <- NULL
  }
  print.default(x)
}


#' Print function for class siftr_matrix 
show.siftr_matrix <- function(sift_mat){
  print(sift_mat)
}

#' Summary for siftr matrix 
#' 
#' @param sift_mat siftr matrix generated using predictFromAlignment
#' 
#' @export
summary.siftr_matrix <- function(sift_mat){
  aln_file =  attr(sift_mat, 'aln_file')
  aln_file_msg = ifelse(is.na(aln_file), ': character vector', 
                        paste( ' alignment file:', aln_file) )
  
  ncol = min(6, ncol(sift_mat))
  cat(sprintf('Siftr matrix object (only showing 5 rows, %s columns)\n', ncol))
  cat(sprintf('generated using%s \n\n', aln_file_msg) )
  
  print(sift_mat[1:6, 1:ncol])
}

#' Generic information content function
#' 
#' @export
infoContent <- function(x){
  UseMethod("infoContent", x)
}

#' Position information content 
#' 
#' @param sift_mat siftr matrix generated using predictFromAlignment
#' 
#' @return numeric vector of information contents
#' @export
infoContent.siftr_matrix <- function(sift_mat){
  attr(sift_mat, 'info_content')
}


#' Generic information content function
#' 
#' @export
diffAA <- function(x){
  UseMethod("diffAA", x)
}

#' Number of unique amino acids per position
#' 
#' @param sift_mat siftr matrix generated using predictFromAlignment
#' 
#' @return integer vector, each value is between 0-20
#' 
#' @export
diffAA.siftr_matrix <- function(sift_mat){
  attr(sift_mat, 'diff_aas')
}


########################################
# Some extra notes <ignore>

#MAXSEQ -> maximum number of sequences allowed <= 400
#RESIDUE_THRESHOLD -> number of amino acids must >= least two 
#THRESHOLD -> sift score threshold 

# # Check we haven't exceeded maximum number of sequences
# # If we have take first N
# if(length(aln_fa) > max_seqs){
#   aln_fa = aln_fa[1:max_seqs]
# }



## RUN PSIBLAST
# $NCBI/blastpgp -d $seq_database -i $query -o $query.out -m 0 -j $iterations -e .0001 -h 0.002 -b 399
# NEWBLAST:
#    psiblast ...


## CONVERT OUTPUT TO FASTA
# $bindir/psiblast_res_to_fasta_dbpairwise $query.out $query.globalX $iterations $query.unfiltered

# query_seq = read_sequence_from_filename (queryfilename);
# 
# read_psiblast_header_until_last (alignmentfp, max_iterations); 
# nseqs = psiblast_pairwise (alignmentfp, alignedseqs,
#                            query_seq->length, option_carve_gaps); \


## Clump output 
# $bindir/clump_output_alignedseq $query.globalX $tmpdir/$pid.clumped .9 0 #>& /dev/null

## Choose sequences via psiblastseedmedian
# $bindir/choose_seqs_via_psiblastseedmedian $query.unfiltered $tmpdir/$pid.clumped $query.selectedclumped 1.0 
# $tmpdir/$pid $median_threshold >& /dev/null

# $bindir/consensus_to_seq $query.selectedclumped $tmpdir/$pid.clumped.consensuskey $tmpdir/$pid.selected >& /dev/null
# 
# 
# $bindir/seqs_from_psiblast_res $query.out $tmpdir/$pid.selected $iterations $query.unfiltered $tmpdir/$pid.alignedfasta $pid
