#' Run BLAST query
#'
#' Perform local BLAST search using requested BLAST type (blastn, blastp., etc.).
#'
#' If `query` or `subject` are the sequnces from [read_fasta()], these are saved into
#' a tempfile and their filenames are passed to the blast.
#'
#' Exactly one of `subject` and `db` must be specified.
#'
#' The output of blast run is stored in a file specified by the `out` argument.
#' The output is then read and returned. If `out = NULL`, the output is stored in a tempfile
#' and will be removed at the end of the current R session.
#'
#' @param query a fasta file or a sequences object.
#' @param subject **optional** fasta file or a sequences object against which to blast,
#' one of subject or db needs to be specified.
#' @param db **optional** a database created from [make_blast_db()]
#' one of subject or db needs to be specified.
#' @param type **optional** type of BLAST to run, one of `blastn`, `blastp`, `blastx`,
#' `tblastn` or `tblastx`.
#' @param out **optional** output path to store BLAST results.
#' @param outfmt **optional** output format, defaults to `0` (Pairwise).
#' @param args **optional** a character vector of additional arguments passed to blast.
#' @param quiet **optional** if TRUE, any blast output to stdout and stderr is suppressed
#' @return results of BLAST run in a requested format.
#'
#' @export
blast = function(
    query,
    subject = NULL,
    db = NULL,
    type = NULL,
    out = NULL,
    outfmt = NULL,
    args = NULL,
    quiet = FALSE
    ){
    write_fasta_temp = function(x){
        file = tempfile()
        write_fasta(x, file)
        file
        }

    types = c("blastn", "blastp", "blastx", "tblastn", "tblastx")
    type = match.arg(type, types)

    if(!xor(is.null(subject), is.null(db)))
        stop("Exactly one of subject or db must be specified.")

    if(is_sequences(query))
        query = write_fasta_temp(query)

    if(is_sequences(subject))
        subject = write_fasta_temp(subject)

    if(is.null(out)){
        out = tempfile()
        }

    args = c(
        "-query", query,
        if(!is.null(db)) c("-db", db),
        if(!is.null(subject)) c("-subject", subject),
        "-out", out,
        if(!is.null(outfmt)) c("-outfmt", outfmt),
        if(!is.null(args)) args
        )

    std = if(isFALSE(quiet)) "" else FALSE
    errcode = system2(type, args, stderr = std, stdin = std)

    if(errcode != 0)
        stop(type, " encountered error.")

    if(!is.null(outfmt) && outfmt %in% c(6,7,10)){
        blast_format(out, outfmt = outfmt)
        } else {
        readLines(out)
        }
    }


#' @rdname blast
#' @export
blastn = function(
    query, subject = NULL, db = NULL,
    out = NULL, outfmt = NULL,
    args = NULL, quiet = FALSE
    ){
    blast(query, subject, db, out, outfmt, args, type = "blastn", quiet = quiet)
    }


#' @rdname blast
#' @export
blastp = function(
    query, subject = NULL, db = NULL,
    out = NULL, outfmt = NULL,
    args = NULL, quiet = FALSE
    ){
    blast(query, subject, db, out, outfmt, args, type= "blastp", quiet = quiet)
    }


#' @rdname blast
#' @export
blastx = function(
    query, subject = NULL, db = NULL,
    out = NULL, outfmt = NULL,
    args = NULL, quiet = FALSE
    ){
    blast(query, subject, db, out, outfmt, args, type= "blastx", quiet = quiet)
    }


#' @rdname blast
#' @export
tblastn = function(
    query, subject = NULL, db = NULL,
    out = NULL, outfmt = NULL,
    args = NULL, quiet = FALSE
    ){
    blast(query, subject, db, out, outfmt, args, type= "tblastn", quiet = quiet)
    }


#' @rdname blast
#' @export
tblastx = function(
    query, subject = NULL, db = NULL,
    out = NULL, outfmt = NULL,
    args = NULL, quiet = FALSE
    ){
    blast(query, subject, db, out, outfmt, args, type= "tblastx", quiet = quiet)
    }


blast_format = function(x, text = NULL, outfmt = NULL){
    # parse only known formats
    if(is.null(outfmt) || !(outfmt %in% c(6, 7, 10))){
        y = if(is.null(text)) readLines(x) else text
        return(y)
        }

    # TODO parse @delim
    # so far we assume the default output
    sep = if(outfmt == 10) "," else "\t"
    x = if(is.null(text)) x else textConnection(text)
    y = utils::read.table(x, sep = sep)
    names(y) = c(
        "query",
        "subject",
        "perc_identity",
        "alignment_length",
        "mismatches",
        "gap_pens",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
        "evalue",
        "bit score"
        )

    y
    }
