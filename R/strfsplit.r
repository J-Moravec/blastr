#' Split string according to first match
#'
#' Split the character vector `x` into two substrings according to the first matched pattern.
#'
#' This is similar to [base::strsplit()], but instead dividing the string along every
#' match of the pattern, the string is split only along the first match.
#' I.e., similar to: `lapply(strsplit(x, split), \(x){list(x[1], paste0(x[-1], collapse = split)`
#' Unlike the above, instead of returning a list of vectors of size two, we return a list
#' of size two, each containing a vector with a length equal to to `x`.
#' This way getting the first or the second element is easier.
#'
#' @param x an input string
#' @param pattern a pattern along which to split
#' @param ... other arguments passed to [base::regexpr()], such as `fixed = TRUE`
#' @return a list of two character vectors, both with the same length as `x`
#'
#' @examples
#' s = strfsplit("foo bar", " ", fixed = TRUE)
#' identical(s, list("foo", "bar"))
#'
#' s = strfsplit("foo", " ", fixed = TRUE)
#' identical(s, list("foo", NA_character_))
#'
#' s = strfsplit(NA, " ", fixed = TRUE)
#' identical(s, list(NA_character_, NA_character_))
#'
#' @export
strfsplit = function(x, pattern, ...){
    match = regexpr(pattern, x, ...)
    match[which(match == -1)] = nchar(x)[which(match == -1)] + 1
    len = attr(match, "match.length")
    len[which(len == -1)] = NA
    list(substring(x, 1, match - 1), substring(x, match + len))
    }
