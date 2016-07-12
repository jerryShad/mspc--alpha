# Bioconductor Package for Chip peak colocalization, parallel intersection

#' @title parse.BedFiles
#' @param bedFiles bed files folder
#' @return GRanges objects
#' 
#' @export
#' @author Julaiti Shayiding

parse.bedFiles <- function(bedFiles) {
  if (length(bedFiles) == 1) {
    if (is.Home_dir(bedFiles)) {
      files <- list.files(path=bedFiles)
      idx <- unlist(sapply(c("bed", "bedGraph"), grep, files))
      idx <- sort(unique(idx))
      files <- files[idx]
      bedFiles <- sub("/$", "", bedFiles)
      res <- paste(bedFiles, files, sep="/")
    }
    else {
      if (!file.exists(bedFiles)) {
        stop("bedFiles folder is not exists in home directory...")
      } 
      else {
        res <- bedFiles
      }
    }
  }
  else {
    if (is.Home_dir(bedFiles[1])) {
      stop("bedFiles should folder with at leat two bed files in it")
    } else {
      res <- bedFiles[file.exists(bedFiles)]
      if (length(res) == 0) {
        stop("bed file does not exists...")
      }
    }
  }
  return(res)
}

# is.Home_dir
#' @description 
#' check whether bed files in home directory

is.Home_dir <- function(dir=.GlobalEnv) {
  if (file.exists(dir) == FALSE)
    return(FALSE)
  return(file.info(dir)$isdir)
}
