# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title readBedAsGRanges
#' @param peakFolder set of bed files to be read
#' @param readAs class type that read peakFile as GRanges objects
#' @return GRanges object
#' @description
#' we encourage user to read sample file as GRanges objects in R
#'
#' @export
#' @importFrom rtracklayer import.bed
#' @author  Julaiti Shayiding

readBedAsGRanges<- function(peakFolder) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- lapply(1:length(files), function(ele_) {
    res <- as(import.bed(files[ele_]), "GRanges")
  })
  return(f.read)
}
