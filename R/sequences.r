#' Read and write sequences in a fasta format
#'
#' Read and write DNA, RNA, or protein sequences in a fasta format.
#' Files compressed with `gzip` (extension `gz`) are supported.
#'
#' No checks are performed, malformed input might lead to malformed output.
#'
#' @param file input or output file with fasta sequences, can be gzipped.
#' @param x a named character vector assumed to be sequences.
#' @param ... arguments passed to the print method (currently unused).
#' @return `read_fasta` returns a named character vector of sequences with class
#' `c("sequences", "character") allowing for type-checking.
#' `write_fasta` doesn't have return value.
#'
#' @examples
#' seq1 = c("A" = "AAA", "B" = "BBB", "C" = "CCC")
#' fasta = tempfile(fileext = ".fasta")
#' writeLines(paste0(">", names(seq1), "\n", seq1), fasta)
#'
#' seq2 = read_fasta(fasta)
#' identical(seq1, seq2)
#'
#' @name fasta
NULL


#' @rdname fasta
#' @export
read_fasta = function(file){
    text = readLines(file)

    starts = grep("^>", text)
    stops = c(starts[-1]-1, length(text))

    sequences = mapply(
        function(start, stop, text){
            seq = text[(start+1):stop]
            seq = gsub("[[:blank:]]*", "", seq)
            paste0(seq, collapse="")
            },
        starts, stops, MoreArgs = list(text)
        )
    names(sequences) = sub("^>", "", text[starts])

    structure(sequences, class = c("sequences", class(sequences)))
    }


#' @rdname fasta
#' @param nchar **optional** a number of characters per line
#' @export
write_fasta = function(x, file = "", nchar = 80){
    if(tools::file_ext(file) == "gz")
        file = gzfile(file)

    if(!is.null(nchar)){
        nchar = as.numeric(nchar)
        x = gsub(
            paste0("(.{1,", nchar, "})"),
            "\\1\n",
            x
            )
        }

    text = paste0(">", names(x), "\n", x)
    # remove the last newline, write lines already adds it
    text[length(text)] = sub("\n$", "", text[length(text)])
    writeLines(text, file)
    }

#' @rdname fasta
#' @export
print.sequences = function(x, ...){
    print(unclass(x))
    }
