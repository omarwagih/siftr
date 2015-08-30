setwd('~/Development/siftr/')

detachPackage <- function(pkg){
  pkg = sprintf("package:%s", pkg)
  res = tryCatch({
    detach(pkg, unload=TRUE, character.only=T, force=T)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(res)
}

.SIFTR_VERSION = '1.0.0'
build_package <- function(){
  # Update description file
  
  require(devtools)
  targz = sprintf('siftr_%s.tar.gz', .SIFTR_VERSION)
  # Move up one directory
  newf = file.path('./build', targz)
  # Compile things
  document('./')
  file.remove(newf)
  system('R CMD BUILD ./')
  re = file.rename(targz, newf)
  # Install package
  system(sprintf('R CMD INSTALL %s', newf))
  # Reload in current environment
  detachPackage('siftr')
}

# Build package
build_package()