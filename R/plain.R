#' Formatting to plain text
#'
#' @param x string: text to be formatted
#' @param ... further arguments passed to or from other methods
#' @export

plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
}
