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
    args = NULL
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

    if(inherits(query, "sequences"))
        query = write_fasta_temp(query)

    if(inherits(subject, "sequences"))
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

    errcode = system2(type, args)

    if(errcode != 0)
        stop(type, " encountered error.")

    readLines(out)
    }


#' @rdname blast
#' @export
blastn = function(query, subject = NULL, db = NULL, out = NULL, outfmt = NULL, args = NULL){
    blast(query, subject, db, out, outfmt, args, type = "blastn")
    }


#' @rdname blast
#' @export
blastp = function(query, subject = NULL, db = NULL, out = NULL, outfmt = NULL, args = NULL){
    blast(query, subject, db, out, outfmt, args, type= "blastp")
    }


#' @rdname blast
#' @export
blastx = function(query, subject = NULL, db = NULL, out = NULL, outfmt = NULL, args = NULL){
    blast(query, subject, db, out, outfmt, args, type= "blastx")
    }


#' @rdname blast
#' @export
tblastn = function(query, subject = NULL, db = NULL, out = NULL, outfmt = NULL, args = NULL){
    blast(query, subject, db, out, outfmt, args, type= "tblastn")
    }


#' @rdname blast
#' @export
tblastx = function(query, subject = NULL, db = NULL, out = NULL, outfmt = NULL, args = NULL){
    blast(query, subject, db, out, outfmt, args, type= "tblastx")
    }
