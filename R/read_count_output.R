#' Read Count Output
#'
#'Generate labelled count matrix
#'
#' @param dir string: path to folder containing files with `name` as prefix
#' @param name string: prefix of the gene/barcode/count files
#' @return Labelled Count matrix object
#' @export

# Function - Generate labelled count matrix
read_count_output <- function(dir, name) {

  dir <- normalizePath(dir, mustWork = TRUE)
  m <- Matrix::readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- methods::as(m, "dgCMatrix")
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}
