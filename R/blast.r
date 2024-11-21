#' Run BLAST query
#'
#' @description
#' Perform local BLAST search using requested BLAST type (blastn, blastp., etc.).
#' Arguments are passed to the NCBI BLAST+ software that must be intalled independently.
#'
#' @details
#' The function `blast` implements the general logic of calling the NCBI BLAST+ using the
#' [system2()] interface. Other functions like `blastn`, `blastp` are calling `blast` with
#' specified `type`. Use them if you want to be more specific in what blast search
#' you are performing.
#'
#' NCBI BLAST+ can process the subject database in a fasta format or in the form of database
#' created using [make_blast_db()]. For this reason, exactly one of `subject` or `db` must
#' be specified. In practice, aside from very small queries, it is almost always more
#' time-efficient to create the database even for one-time queries.
#
#' If `query` or `subject` are the sequnces, such as from [read_fasta()], these are saved into
#' a tempfile and their filenames are passed to NCBI BLAST+.
#'
#' The output of blast run is stored in a file specified by the `out` argument.
#' The output is then read and returned. If `out = NULL`, the output is stored in a tempfile
#' and will be removed at the end of the current R session.
#'
#' If the `outfmt` is one of the tabular formats (6, 7, or 10), the NCBI BLAST+  is parsed
#' and converted into a `data.frame`. As of now, field customization is not supported.

#' @param query A fasta file or a sequences object.
#' @param subject Fasta file or a sequences object against which to blast,
#' one of `subject or `db` needs to be specified.
#' @param db A database created from [make_blast_db()],
#' one of `subject` or `db` needs to be specified.
#' @param type Type of BLAST to run, one of `blastn`, `blastp`, `blastx`,
#' `tblastn` or `tblastx`, defaults to `blastn`.
#' @param out An output path to store BLAST results, the default is to not store any results.
#' @param outfmt BLAST output format, defaults to `0` (Pairwise).
#' @param args AA character vector of additional arguments passed to [blast()].
#' @param quiet If `TRUE`, any blast output to stdout and stderr is suppressed,
#' use this when using [blast()] within a pipeline and you want to reduce the text output
#' to console. This is also make potential error harder to debug. The default is `FALSE`.
#' @return a character vector or a data.frame (for the tabular formats `outfmt = 6`, `7`, and
#' `10`) containin the results of BLAST run.
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
    types = c("blastn", "blastp", "blastx", "tblastn", "tblastx")
    type = match.arg(type, types)

    if(!xor(is.null(subject), is.null(db)))
        stop("Exactly one of subject or db must be specified.")

    if(is_sequences(query))
        query = write_fasta_temp(query)

    if(is_sequences(subject))
        subject = write_fasta_temp(subject)

    if(identical(tools::file_ext(query), "gz")){
        y = tools::file_path_sans_ext(query) |> basename()
        query = gunzip(query, file.path(tempdir(), y), keep = TRUE)
        }

    if(identical(tools::file_ext(subject), "gz")){
        y = tools::file_path_sans_ext(subject) |> basename()
        subject = gunzip(subject, file.path(tempdir(), y), keep = TRUE)
        }

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
    errcode = system2(type, args, stderr = std, stdout = std)

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
