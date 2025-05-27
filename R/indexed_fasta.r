#' Read indexed fasta
#'
#' Read a fasta file by using pre-computed index `.fai` using random access.
#'
#' The index `.fai` is a precomputed file containing starts and lengths of sequences within
#' the fasta file. This allows obtaining specific sequences from the fasta file without
#' having to first read the whole file sequentially.
#'
#' `read_indexed_fasta()` reads the index file and creates a `indexed_fasta` object.
#'
#' The `indexed_fasta` object is internally represented as a list with a `file` the input
#' fasta and the index `index`. Functions such as `[`, `[[` and `names()` are defined
#' to make working with this object more pleasant and intutitive.
#' To view the internal structure of `obj` of class `indexed_fasta`, use `unclass(obj)`.
#'
#' `[` is used to subset the `indexed_fasta` by subsetting the index.
#'
#' `[[` is used to read one or more sequences from `indexed_fasta` using the [base::seek()].
#' It returns object of the class `sequences` (see `blastr::sequences`)
#'
#' Methods `names()`, `names<-()` allow getting or setting names,
#' `length()` allows for obtaining the number of sequences within the `indexed_fasta` object.
#'
#' @param x a fasta file
#' @param index an index to the corresponding fasta file
#' @param value in `[`, `[[`, or `$`, either a numeric index or a character string
#' which will be subsetted
#' @param ... other arguments passed to or from functions
#' @return an object of class `indexed_fasta`
#'
#' @name indexed_fasta
NULL

# internal
read_fasta_index = function(x){
    y = utils::read.table(x, sep = "\t", header = FALSE)
    names(y) = c("name", "length", "offset", "linebases", "linewidth")

    # calculate number of lines to read with readLines
    y$lines = ceiling(y$length / y$linebases)
    y
    }


#' @rdname indexed_fasta
#' @export
read_indexed_fasta = function(x, index = NULL){
    if(is.null(index))
        index = paste0(x, ".fai")
    if(!file.exists(index)){
        # TODO try to create it
        stop("File \"", index, "\" not found.")
        }
    if(is_gzip(x))
        stop("gzip doesn't support random access.")
    index = read_fasta_index(index)
    structure(list("index" = index, file = x), class = "indexed_fasta")
    }


#' @rdname indexed_fasta
#' @export
"[.indexed_fasta" = function(x, value){
    cl = class(x)
    x = unclass(x)
    if(is.character(value))
        value = match(value, x$index$name)
    x$index = x$index[value,]
    structure(x, class = cl)
    }


#' @rdname indexed_fasta
#' @export
"[[.indexed_fasta" = function(x, value){
    x = unclass(x)

    if(is.character(value))
        value = match(value, x$index$name)

    x$index = x$index[value,]

    con = file(x$file, "r")
    on.exit(close(con), add = TRUE)

    n = length(value)
    seqs = vector("character", n)
    for(i in seq_len(n)){
        seek(con, x$index$offset[i])
        seqs[i] = paste0(readLines(con, x$index$lines), collapse = "")
        }
    names(seqs) = x$index$name
    structure(seqs, class = c("sequences", "character"))
    }


#' @rdname indexed_fasta
#' @export
print.indexed_fasta = function(x, ...){
    x = unclass(x)
    cat(sprintf(
        "indexed fasta: %d sequences of average size %.2f\n",
        nrow(x$index), mean(x$index$length)), sep = "")
    }


#' @rdname indexed_fasta
#' @export
length.indexed_fasta = function(x){
    nrow(unclass(x)$index)
    }


#' @rdname indexed_fasta
#' @export
names.indexed_fasta = function(x){
    unclass(x)$index$name
    }


#' @rdname indexed_fasta
#' @export
"names<-.indexed_fasta" = function(x, value){
    cl = class(x)
    x = unclass(x)
    x$index$name = value
    structure(x, class = cl)
    }

#' @export
"$.indexed_fasta" = function(x, value){
    if(value %in% names(x)){
        x[[value]]
        } else {
        NULL
        }
    }

# Overrwide some setters

#' @export
"$<-.indexed_fasta" = function(x, index, value){
    stop("\"$<-\" is undefined for \"indexed_fasta\"")
    }

#' @export
"[<-.indexed_fasta" = function(x, index, value){
    stop("\"[<-\" is undefined for \"indexed_fasta\"")
    }

#' @export
"[[<-.indexed_fasta" = function(x, index, value){
    stop("\"[[<-\" is undefined for \"indexed_fasta\"")
    }
