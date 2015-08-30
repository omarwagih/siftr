# For help see <http://sift.bii.a-star.edu.sg/www/SIFT_help.html>

#' Remove identical sequences in a list alignment 
#' 
#' Given an alignment of sequences (as a list), 
#' Remove sequences that have 100% identity to query (first element)
#' 
#' @param aln alignment list of sequences
#' @return alignment list with identical sequences removed
removeIdenticalSeqs <- function(aln){
  id = sapply(aln, percentIdentity, aln[[1]]) >= 100 #%
  # Always keep query
  id[1] = FALSE
  return(aln[!id])
}

#' Compute percent identity between two sequences
#' 
#' @param comp sequence being compared
#' @param seq suery sequence
#' 
#' @return percent identity between both sequences
percentIdentity <- function(comp, seq){
  ident = sum(seq == comp)
  ident = ident - sum(comp == 'X')
  round( ident / length(seq ) * 100)
}

#' SIFT version of median
#' 
#' The way medians are computed in sift C code is a bit different
#' This implements the median function used by sift
#' 
#' @param arr vector of numeric values
#' 
#' @return median value
medianOfArray <- function(arr){
  n = length(arr)
  
  if(n == 1) return(array)
  
  arr = sort(arr)
  z = length(arr)/2
  
  if( (n %% 2) == 0 ){
    return(arr[z])
  } else {
    return((arr[floor(z)] + arr[ceiling(z)])/2)
  }
}

#' Check lengths of sequences are equal
#' 
#' This checks that the lengths of each sequence in the 
#' alignment is equal to the length of the query
#' 
#' @param aln alignment list of sequences
#' 
#' @return will throw error if there is at least one unequal sequence length
#' otherwise, will return length of the query
checkUnequalSequenceLengths <- function(aln){
  # Get lengths of all sequences 
  seq_len = sapply(aln, length)
  # Stop if they're not all the same as the query length
  if(! all(seq_len == seq_len[1])) stop('Unequal length of sequences')
  
  # Return number of positions
  seq_len[1]
}

#' Compute matrix of letters for the alignment
#' 
#' This takes all n sequences, of m positions and constructs
#' a matrix containing sequence amino acid x position
#' amino acids are represented as integers 1-20, 21 non-AA character
#' 
#' @param aln alignment list of sequences
#' @param alphabet character vector of alphabet (amino acids used)
#' @param npos number of positions in the alignment
#' 
#' @return matrix of amino acid indices for each sequence/position
letterMatrix <- function(aln, alphabet, npos){
  # Unlist all residues
  split = unlist( aln, F, F )
  
  # Convert to numeric
  split = as.numeric( factor(split, levels = alphabet) )
  split[is.na(split)] = length(alphabet) + 1
  
  # Reformat into a matrix with correct number of rows/columns
  aln_mat = t( matrix(split, npos, length(split)/npos ) )
  
  aln_mat
}

#' Get median information content for each position
#' 
#' @param aln_mat alignment matrix generated from letterMatrix
#' @param frq precomputed amino acid frequency in the proteome
#' @param cores number of cores for parallel processing (default: 1)
#' 
#' @return Numeric vector with median information content for each position
medianInfoContent <- function(aln_mat, frq, cores=1){
  # Used for information content calculation
  ln2 = logb(2.0)
  hmax = logb(20) #20 amino caids
  
  # Position profiles / what sequences to keep for each position
  position_profiles = apply(aln_mat, 2, function(x) (x<21)+0)
  
  # Convert to sequence
  position_profiles_seq = apply(position_profiles, 2, paste0, collapse='')
  
  # Get unique profiles / avoid running duplicates
  unique_profiles = which(!duplicated(position_profiles_seq, MARGIN = 2))
  
  # Match all positions back to the unique one
  info_content_ind = match(position_profiles_seq, 
                           position_profiles_seq[unique_profiles])
  
  if(cores == 1){
    # Single core using vapply (faster than sapply)
    median_info_content = vapply(unique_profiles, FUN.VALUE = double(1), 
                                 positionInfoContent, aln_mat, frq, FALSE)
  }else{
    # Parallel process using mclapply
    median_info_content = mclapply(unique_profiles, mc.cores = cores, 
                                   FUN = positionInfoContent, aln_mat, frq, FALSE)
    median_info_content = simplify2array(median_info_content)
  }
  
  #positionInfoContent(1, aln_mat, F)
  median_info_content = median_info_content[info_content_ind]
}

#' Get median information content for a given position
#' 
#' Get sequences that have no gap or X in position i, 
#' compute median information content accross all positions
#' 
#' @param position alignment position of interest 
#' @param aln_mat alignment matrix generated from letterMatrix
#' @param frq precomputed amino acid frequency in the proteome
#' @param small_sample_correction TRUE if small sample correction is applied 
#' to information content calculation (not implemented currently, sift has it but doesn't use it)
#' @return information content value between 0 (conserved) to log2(20) (least conserved) 
positionInfoContent <- function(position=1, aln_mat, frq, small_sample_correction=F){
  # Used for information content calculation
  ln2 = logb(2.0)
  hmax = logb(20) #20 amino caids
  
  # Number of positions
  npos = ncol( aln_mat )
  
  # Which sequences have valid residues
  ind = aln_mat[,position] < 21
  
  # Number of sequences we're computing info content on
  nseqs = sum(ind)
  
  # Handle case where we have 0 amino acids (gaps or non amino acids)
  if(nseqs == 0) return(NA)
  # Handle case where we have one sequence
  if(nseqs == 1) return(log2(20))
  
  # Extract sub matrix
  aln_mat_curr = aln_mat[ind,]
  
  # writeLines(sprintf('position = %s', position))
  
  # Get pb weights
  ptm <- proc.time()
  pb_weights = pbWeights(aln_mat_curr)
  # writeLines( sprintf('pbWeights took %f seconds', (proc.time() - ptm)['elapsed']) )
  
  
  # Make PWM
  ptm <- proc.time()
  # Make PWM for sequences with information
  pwm = makePWM( aln_mat_curr, pb_weights, frq )
  # writeLines( sprintf('makePWM took %f seconds', (proc.time() - ptm)['elapsed']) )
  
  ptm <- proc.time()
  
  # Compute information content for each position
  info_per_pos = sapply(1:npos, function(i){
    pos_data = pwm[,i]
    dtemp = pos_data / sum(pos_data, na.rm=T)
    
    r = hmax
    r = r + sum( dtemp * logb(dtemp), na.rm=T )
    
    r / ln2
  })
  
  
  # writeLines( sprintf('infoContent took %f seconds', (proc.time() - ptm)['elapsed']) )
  
  m = medianOfArray(info_per_pos)
#   writeLines(sprintf('median_info=%s', m))
#   writeLines('_______')
  m
}


#' Compute a weight for each amino acid at each position
#' 
#' These weights are used to compute a sequence weight 
#' denoting how much information each sequence is contributing
#' 
#' @param aln_mat alignment matrix generated from letterMatrix
#' @param return_all if TRUE, returns intermediate variables computed (pfm, diff_aas, etc), otherwise just weight matrix
#' 
#' @return numeric matrix of n sequences x m positions, with a weight for each cell
pbWeights <- function(aln_mat, return_all=FALSE){
  # Number of sequences/positions
  npos = ncol( aln_mat )
  nseq = nrow( aln_mat )
  
  # Construct data table with position, amino acid, and freq (to be grouped)
  pb_data = data.table(pos=rep(1:npos, each = nseq), 
                       aa=as.integer(aln_mat), freq=1)
  
  # Compute amino acid frequency (group by pos and aa)
  pb_counts = pb_data[, lapply(.SD, sum), by=c('pos', 'aa')]
  
  # Pick valid amino acids non gaps and non X (1 to 20)
  pb_counts_valid = pb_counts[pb_counts$aa < 21,]
  # Compute number of different amino acids per position
  diff_aas = ftable(factor(pb_counts_valid$pos, 1:npos))
  
  # Compute pb_weight
  pb_counts$pb = 1 / (pb_counts$freq * diff_aas[ pb_counts$pos ])
  pb_counts = as.matrix( pb_counts[pb_counts$aa < 21,] )
  
  # Construct pfm / pb matrix from data and assign counts
  pb_mat = pfm = matrix(0, 21, npos)
  pfm[ pb_counts[,2:1] ] = pb_counts[,3]
  pb_mat[ pb_counts[,2:1] ] = pb_counts[,4]
  
  # Loop over each position and get pb_weights per sequence (used to compute sequence weight)
  dat = vapply(1:ncol(aln_mat), FUN.VALUE = numeric( nseq ), FUN=function(pos){
    pb_mat[aln_mat[,pos],pos] # pb weight
  })
  
  
  if(return_all){
    return( list(pb_mat = dat, pfm = pfm, 
                 diff_aas = as.integer(diff_aas), 
                 pb_counts=pb_counts, npos=npos, nseq=nseq) )
  }
  
  dat
}



#' Make PWM using pb weights and amino acid frequencies
#' 
#' Computes a weight matrix of 20 amino acids x n positions
#' Each cell is the frequency of the corresponding amino acid in that position * weight
#' 
#' @param aln_mat alignment matrix generated from letterMatrix
#' @param pb_weights weights computed from pbWeights function
#' @param frq precomputed amino acid frequency in the proteome
#' 
#' @return weight matrix of 20 amino acids x number of positions 
makePWM <- function(aln_mat, pb_weights, frq){
  # Number of positions and sequences
  npos = ncol( aln_mat )
  nseq = nrow( aln_mat )
  
  # Get weights for each sequence
  seq_weights = rowSums(pb_weights, na.rm=T)
  
  # Data for building pwm: position, amino acid and sequence weight (from pb_weights)
  pwm_df = data.frame(pos=rep(1:npos, each = nseq),
                  aa=as.integer(aln_mat), 
                  seq_w = rep(seq_weights, times = npos))
  
  # Assign background amino acid frequencies
  pwm_df$bg_freq = frq[ pwm_df$aa ]
  
  # Compute new weight
  pwm_df$new_weight = pwm_df$seq_w / pwm_df$bg_freq
  
  # Construct data table (much faster to aggregate)
  dt = as.data.table( pwm_df[,c('pos', 'aa', 'new_weight')] )
  
  # Aggregate weights table by position and amino acid
  dt = dt[, lapply(.SD, sum), by=c('pos', 'aa')]
  
  # Choose only valid amino acids
  dt = as.matrix(dt[dt$aa < 21,])
  
  # Construct pwm and assign weights
  pwm = matrix(0, 21, npos)
  pwm[ dt[,2:1] ] = dt[,3]
  pwm
}



#' Normalize pb weights matrix
#' 
#' This takes wight matrix produced by pbWeights and normalizes 
#' by number of sequences in alignment / sum of all weights
#' 
#' @param weights weights computed from pbWeights function
#' @return normalized pb weight matrix
normalizeWeights <- function(weights){
  w_sum = sum(weights)
  num_sequences = nrow(weights)
  factor = num_sequences / w_sum
  weights = factor * weights
  return(weights)
}

#' Normalizes weights computed by pbWeights
#' 
#' This normalizes the pbWeights weight matrix using normalizeWeights 
#' then groups by amino acid for 20 x number of position matrix
#' This is used to compute sume of weights for each position
#' 
#' @param pb_weights weights computed from pbWeights function
#' @param aln_mat alignment matrix generated from letterMatrix
#' @return list of two elements: (1) matrix of normalized weights 
#' grouped by amino acids (20 x number of positions)
#' (2) numeric vector of the sum of weights per positions
pbWeightsNorm <- function(pb_weights, aln_mat){
  # Number of positions/sequences
  npos = ncol(aln_mat)
  nseq = nrow(aln_mat)
  
  # Normalize sequence weights
  pb_weights_norm = normalizeWeights(pb_weights$pb_mat)
  
  # Get new sequence weights
  seq_weights = rowSums(pb_weights_norm, na.rm=T)
  
  # Construct data table (for fast aggregation)
  dt = data.table(pos = rep( 1:npos, each = nseq ), 
                  aa = as.integer( aln_mat ), 
                  freq = 1,
                  seq_w = rep( seq_weights, times = npos ))
  
  # Group sum weights by position and amino acid
  dt = dt[, lapply(.SD, sum), by=c('pos', 'aa')]
  
  # Remove non amino acid residues
  dt_mat = as.matrix( dt[dt$aa < 21] )
  
  # Make into weight matrix
  # cnt (list) is now sum_seq_weights (matrix, columns) 
  grouped_seq_weights = matrix(0, 21, npos)
  grouped_seq_weights[ dt_mat[,2:1] ] = dt_mat[,4]
  
  # totcnt -> sum_seq_weights
  sum_seq_weights = apply(grouped_seq_weights, 2, sum)
  
  list(grouped = grouped_seq_weights[1:20,], 
       sum = sum_seq_weights)
}

#' Compute epsilon values 
#' 
#' Used for computing final sift scores
#' 
#' @param pb_weights see pbWeights
#' @param pb_weights_norm see pbWeightsNorm
#' @param rank_matrix precomputed rank matrix for amino acids (BLOSUM62)
getEpsilonValues <- function(pb_weights, pb_weights_norm, rank_matrix){
  # similarity_dependent_scale_0
  dtemp2 = apply(pb_weights_norm$grouped, 2, function(pos_data){
    # Get heaviest weight amino acid
    max_aa = which.max(pos_data)
    # Get rank with other amino acids
    rank = rank_matrix[max_aa,]
    w = pos_data / sum(pos_data)
    sum( rank * pos_data / sum(pos_data) )
  })
  
  # Compute epsilon values
  # different amino acids per position
  diff_aas = pb_weights$diff_aas
  eps = integer(length(diff_aas))
  
  ind = diff_aas > 1
  eps[ind] = exp(dtemp2[ind])
  eps
}


#' Computes final matrix of sift scores returned to the user
#' 
#' @param pb_weights see pbWeights
#' @param pb_weights_norm see pbWeightsNorm
#' @param eps epsilon values
#' @param diri precomputed dirichelet components, used to compute pseudocounts
#' 
#' @return numeric matrix with amino acid rows x number of positions containing sift scores
makeSiftMatrix <- function(pb_weights, pb_weights_norm, eps, diri){
  # Total number of amino acids per columns (totreg)
  num_aa = apply(pb_weights$pfm, 2, sum)
  npos = pb_weights$npos
  
  # Seq weights
  sum_seq_weights = pb_weights_norm$sum
  grouped_seq_weights = pb_weights_norm$grouped
  
  # Construct final sift matrix of siftr scores
  sift_matrix = vapply(1:npos, FUN.VALUE = double(20), FUN = function(i){
    d = pseudo_diri(pos=i, diri, grouped_seq_weights, sum_seq_weights)
    
    fi = grouped_seq_weights[,i] + (eps[i] * d$pseudocounts)
    ind = sum_seq_weights[i] > d$totreg
    if(any(ind)) fi = fi / (sum_seq_weights[i] + eps[i])
    
    fi / max(fi)
  })
  
  sift_matrix
}

#' Helper function for pseudo_diri, adds two logs
#' 
#' @param lx first value
#' @param ly second value
addLogs <- function(lx, ly) {
  if (lx > ly) return (lx + log(1.0 + exp(ly - lx)))
  else return (ly + log(1.0 + exp(lx - ly)));
}

#' Get dirichelet pseudocounts for position i
#' 
#' @param pos position
#' @param diri precomputed dirichelet components, used to compute pseudocounts
#' @param grouped_seq_weights see pbWeightsNorm
#' @param sum_seq_weights see pbWeightsNorm
#' 
#' @return pseudocounts for amino acids of position i
pseudo_diri <- function(pos, diri, grouped_seq_weights, sum_seq_weights){
  
  # Get diri objects
  diri_alpha = diri$alpha
  altot = diri$altot
  q = diri$q
  
  # Number of components
  ncomp = length(altot)
  
  # Weights of amino acids for this position (cnt)
  seq_weights_pos = grouped_seq_weights[1:20,pos]
  # Sum of seq_weights_pos (totcnt)
  sum_seq_weights_pos = sum_seq_weights[pos]
  
  # Keep weights that have a weight
  nonzero_weights = seq_weights_pos > 0
  seq_weights_pos = seq_weights_pos[nonzero_weights]
  
  # compute equation (3), Prob(n|j)
  probn = sapply(1:ncomp, function(j){
    
    lg = lgamma( sum_seq_weights[pos] + 1.0 ) + lgamma(altot[j])
    lg = lg - lgamma(sum_seq_weights[pos] + altot[j])
    
    dtemp = lgamma(seq_weights_pos + diri_alpha[nonzero_weights, j])
    dtemp = dtemp - lgamma(seq_weights_pos + 1.0)
    dtemp = dtemp - lgamma(diri_alpha[nonzero_weights, j])
    
    lg + sum(dtemp)
  })
  
  denom = log(q[1]) + probn[1]
  dtemp = log(q[-1]) + probn[-1]
  for(j in 1:(ncomp-1)) denom = addLogs(denom, dtemp[j])
  
  probj = log(q) + probn - denom
  
  aa_reg = sapply(1:20, function(aai){
    reg = exp(probj) * diri_alpha[aai,]
    sum(reg)
  })
  aa_totreg = sum(aa_reg)
  
  aa_reg = aa_reg / aa_totreg
  pseudocounts = aa_reg
  
  # Return pseudo counts
  return(list(pseudocounts=pseudocounts, 
              totreg=aa_totreg))
}


