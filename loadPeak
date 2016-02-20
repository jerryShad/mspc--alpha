# Bioconductor package for parsing Chip peak over multiple sample

#' @title loadPeak
#' 
#' @param peakFile input bed file to load in R session
#' @param verbose logical indicator that currently non return
#' @return GRanges object with all genomic ranges with required metadata column
#' @export
#' @importFrom rtracklayer import.bed
#' @author Julaiti Shayiding


loadPeak <- function(peakFile, verbose=FALSE){
  if(is(peakFile, "GRanges")){
    peak.gr <- peakFile
  }
  else if(file.exists(peakFile)){
    if(verbose){
      cat(">> loading bed file...\t\t\t",
          format(Sys.time(), "%DD-%MM-%YY %X"), "\n")
      peak.gr <- readSample(peakFile, readAs="GRanges")
    }
  }
  else{
    stop("peakFile should be desired GRanges object")
  }
  return(peak.gr)
}


## readSample
readSample<- function(peakFile, readAs="GRanges") {
  readAs <- match.arg(readAs, "GRanges")
  if(readAs=="GRanges"){
    peakFile <- import.bed(peakFile)
    peak.gr <- as(peakFile, "GRanges")
  }
  return(peak.gr)
}
