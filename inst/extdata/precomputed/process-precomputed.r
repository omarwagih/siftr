
# SIFT Load frequencies
load_frequencies <- function(file='default.amino.frq'){
  li = readLines(file)
  li = li[grepl('[A-Z]\\)$', li)]
  frequency = sapply(strsplit(li, '\\s+'), function(i) as.numeric(i[1]))
  names(frequency) = sapply(strsplit(li, '\\s+|\\)'), function(i) i[length(i)])
  return(frequency)
}

load_diri_file <- function(file='default.diri'){
  # This is the order of the data in the diri file, do not change
  AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  
  li = readLines(file)
  # Find different mixtures 
  nmix = li[grepl('Mixture= ', li)]
  nmix = as.numeric( sapply(strsplit(nmix, '\\s+'), function(i) i[2]) )
  
  # And their alpha values
  alpha = li[grepl('Alpha= ', li)]
  alpha = sapply( strsplit(alpha, '\\s+'), function(i) as.numeric(i[-1])) 
  altot = alpha[1,]
  alpha = alpha[-1,]
  rownames(alpha) = AA
  
  alpha_norm = t(apply(alpha, 1, function(a) a / altot))
  bg_freq = sapply(1:nrow(alpha), function(i) sum(nmix * alpha[i,] / altot) )
  freq_to_bg = sapply(1:ncol(alpha), function(i){
    alpha[,i] / (altot[i] * bg_freq)
  })
  
  return(list(freq_to_bg=freq_to_bg, 
              bg_freq=bg_freq, 
              alpha_norm=alpha_norm, 
              alpha=alpha, 
              altot=altot,
              q=nmix))
}

load_qij <- function(file='default.sij'){
  li = readLines(file)
  
  # Remove comments and split each line
  li = li[!grepl('^#', li)]
  li = sapply(strsplit(li, '\\s+'), function(i) i[i != ''])
  
  # Construct header
  hd = li[[1]]
  n = length(hd)
  vals = t(sapply(li[-1], function(line){
    z = as.numeric(line)
    c(z , rep(NA, n-length(z)))
  }))
  
  # Get NAs and reflect
  na = which(is.na(vals), arr.ind=T)
  vals[na] = vals[na[,2:1]]
  
  # Compute the marginal probabilities
  marg = apply(vals, 1, sum)
  total = unname( sum(marg) )
  rownames(vals) = hd
  colnames(vals) = hd
  
  return(list(vals=vals, 
              marg=marg,
              total=total))
}
NAA=20
construct_rank_matrix <- function(){
  rank_Qij = load_qij()
  vals = rank_Qij$vals
  
  vals = rbind(rep(-1, nrow(vals)), vals)
  aalist = lapply(1:nrow(vals), function(i) c(c(GAP=-50), vals[i,]) )
  
  # Make rank matrix
  rank_matrix = matrix(0, NAA+1, NAA+1)
  for(aa_original in 1:length(aalist)){
    srt = sort( aalist[[aa_original]], decreasing=T)
    sorted_aa = names(srt)
    mch = match(sorted_aa, colnames(vals)) 
    mch[is.na(mch)] = 21
    rank_matrix[aa_original,][mch] = 1:(NAA+1)
  }
  rank_matrix = cbind(rank_matrix[,21],rank_matrix[,-21])
  rownames(rank_matrix) = c('GAP', colnames(vals))
  colnames(rank_matrix) = rownames(rank_matrix)
  return(rank_matrix)
}

# Get into directory
wd <- setwd('~/Development/siftr/precomputed/')

frq = load_frequencies('default.amino.frq')
frq['OTHER'] = NA
saveRDS(frq , 'aa-freq.rds')

rnk = construct_rank_matrix()[-1,-1]
aa_ordered = sort( rownames(rnk) )
rnk = rnk[aa_ordered, aa_ordered]
saveRDS(rnk , 'rank-matrix.rds')

diri = load_diri_file()
# Correct order of amino acids
aa_ordered = order( rownames(diri$alpha) )
diri$alpha = diri$alpha[aa_ordered,]
diri$alpha_norm = diri$alpha_norm[aa_ordered,]
diri$freq_to_bg = diri$freq_to_bg[aa_ordered,]
diri$bg_freq = diri$bg_freq[aa_ordered]
saveRDS(diri , 'diri.rds')

# Set original directory
setwd(wd)