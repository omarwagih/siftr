<img src="https://cdn.rawgit.com/omarwagih/siftr/master/inst/extdata/images/siftr-logo-flat.svg" alt="rmimp logo" width="250px"><br>
<h3>Predicting the impact of mutations on protein function</h3>
===============================================================

## Installation

To install siftr, first make sure your R version is at least R 3.0. You can check by typing the following into your R console:

```r
R.Version()$major
```

Install and load `devtools` package:

```r
install.packages("devtools")
library("devtools")
```

Download and install the `siftr` package from github:

```r
install_github("omarwagih/siftr")
```

Load the `siftr` package

```r
library("siftr")
```

## Running siftr on sample data:

```r
# Get the path to the sample amino acid alignment
# Alternatively you can pass in the alignment as a vector
# e.g. c('SSSS', 'STTT', 'SYSY', 'ASST', 'KKHS')
sample_fa = system.file("extdata", "P39709.alignedfasta", package = "siftr")


# Compute sift scores
# You can adjust the number of cores for faster processing
sift_mat = predictFromAlignment(sample_fa, cores=1)

# Summary of the sift matrix
summary(sift_mat)
```

The generated matrix contains the sift scores. Each row is the substituted amino acid and each column is the position in the alignment. 

Then, filter the matrix for significant scores (< 0.05), and output in a table format:

```r
# Keep predictions with sift scores < 0.05 
filt = filterPredictions(sift_mat, score_thresh = 0.05)

# Display resulting table
print( head(filt) )
```

That's all there is to it!

For further documentation see `?predictFromAlignment` and `?filterPredictions`


## Contact
If you have any feedback, suggestions or questions, please drop me a line at (wagih(at)ebi.ac.uk) or open an issue on github.