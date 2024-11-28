#' Match names
#'
#' Find an elements of `x` in the first space-separated element of `where`.
#'
#' The typical use of this function is to match sequence names.
#' In fasta files, sequence names often consist of the sequence name itself and an additional information about
#' the sequence.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' match_names("foo", "foo bar")
#' match_names(c("foo", "bar"), c("bar", "foo"))
#' match_names(c("foo", "bar") # throws error: foo not found
#' match_names("foo", c("foo", "foo")) # throws error: foo found twice
#' }
match_names = function(x, where){
    where = strfsplit(where, " ", fixed = TRUE)[[1]]
    y = lapply(x, \(y) which(y == where))
    names(y) = x

    # sanity check
    l = lengths(y)
    if(any(l < 1))
        stop("Some sequence names were not found: ", paste0(x[l < 1], collapse = "; "))
    if(any(l > 1))
        stop("Non unique match for sequence names: ", paste0(x[l > 1], collapse = "; "))

    unlist(y)
    }
