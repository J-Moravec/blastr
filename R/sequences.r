#' Read and write sequences in a fasta format
#'
#' Read and write DNA, RNA, or protein sequences in a fasta format.
#' Files compressed with `gzip` (extension `gz`) are supported.
#'
#' No checks are performed, malformed input might lead to malformed output.
#'
#' `read_fasta()` and `write_fasta()` read and write sequences in a fasta format.
#'
#' `write_fasta_temp()` writes a fasta to a `tempfile()` and returns to file.
#' Used internally to convert sequences into a file for NCBI+ BLAST calls.
#'
#' `is_sequences()` tests if the passed object are sequences. In R, operations
#' tend to typically not preserve attributes, like class. `in_sequences()` thus
#' rans several tests instead, but these are only heuristics. If you want to be sure,
#' reassign the class attribute such as with `attr(x, "class") = "sequences"`.
#'
#' @param file input or output file with fasta sequences, can be gzipped.
#' @param x a named character vector assumed to be sequences.
#' @param ... arguments passed to the print method (currently unused).
#' @return
#' `read_fasta` returns a named character vector of sequences with class
#' `c("sequences", "character") allowing for type-checking. Note that attributes
#' like `class` are typically not preserved during various operations (e.g., `head()`).
#' `write_fasta` doesn't have return value and is run for its side effects.
#' `is_sequences()` returns a logical value.
#' `write_fasta_temp()` returns a path to a tempfile.
#'
#' @examples
#' seq1 = c("A" = "AAA", "B" = "BBB", "C" = "CCC")
#' fasta = tempfile(fileext = ".fasta")
#' writeLines(paste0(">", names(seq1), "\n", seq1), fasta)
#'
#' seq2 = read_fasta(fasta)
#' identical(seq1, seq2)
#'
#' @name sequences
#' @aliases fasta
NULL


# For the purpose of printing, or writing, it is faster
# to just insert newline instead of spliting text into
# list of nchar long, see commented out alternative
# implementation
wrap = function(x, nchar){
    gsub(
        paste0("(.{1,", nchar, "})"),
        "\\1\n",
        x
        )
    }

#wrap = function(x, nchar = 80){
#    lapply(x, \(s){
#        i = seq(1, nchar(s), by = nchar)
#        s = substring(s, i, i + nchar - 1)
#        s
#        })
#    }


#' @rdname sequences
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


#' @rdname sequences
#' @param nchar **optional** a number of characters per line
#' @export
write_fasta = function(x, file = "", nchar = 80){
    if(tools::file_ext(file) == "gz")
        file = gzfile(file)

    if(!is.null(nchar)){
        nchar = as.numeric(nchar)
        x = wrap(x, nchar)
        }

    text = paste0(">", names(x), "\n", x)
    # remove the last newline, write lines already adds it
    text[length(text)] = sub("\n$", "", text[length(text)])
    writeLines(text, file)
    }


#' @rdname sequences
#' @export
is_sequences = function(x){
    # NULL is not sequences (or path)
    if(is.null(x)) return(FALSE)

    # many operations remove attributes (such as class)
    # e.g.,: head(x) will lose the class of x
    if(inherits(x, "sequences")) return(TRUE)

    y = utils::head(x)
    # sequences don't have "/", "\" or "\\ paths separators, paths do
    if(any(grepl("/|\\|\\\\", y))) return(FALSE)

    # paths are typically shorter than sequences
    # length of path is essentially unlimited (OS and filesystem specific)
    # but length of file-names is typically shorter than 256 (OS and filesystem specific)
    # since we don't have path separator (previous test), we test for filename length
    if(any(nchar(y) > 256)) return(TRUE)

    # sequences don't have extension
    if(any(tools::file_ext(y) != "")) return(FALSE)

    TRUE
    }


#' @rdname sequences
#' @export
print.sequences = function(x, nchar = 70, ...){
    if(is.numeric(nchar)){
        cat(paste0("[", seq_along(x), "] ", names(x), "\n", wrap(x, nchar)), sep = "\n")
        } else {
        print(unclass(x))
        }
    }


#' @export
"[.sequences" = function(x, i){
    structure(unclass(x)[i], class = class(x))
    }

#' @rdname sequences
#' @export
write_fasta_temp = function(x){
    file = tempfile()
    write_fasta(x, file)
    file
    }
