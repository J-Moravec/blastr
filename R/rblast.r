#' Reciprocal BLAST
#'
#' Perform BLAST reciprocal best hit (RBH) search.
#'
#' BLAST reciprocal best hit search will first perform a forward search, take the best hits,
#' and then performs a backward search using these best hits. If the sequence used for
#' the forward search are identified as a best hit in the backward search, the matched
#' subject are reported.
#'
#' @param x id or names of the sequences to be blasted, id is the first space-separated part
#' of the sequence name and is intepreted by blast as query name
#' @param query,subject sequences of the query and target organisms
#' @param type a type of blast, either nucleotide `blastn` or protein `blastp` blast
#' @param query_db,subject_db **optional** pre-existing query and subject blast dabatases
#' @param keep **optional** keep blast searches, these are returned together with normal output
#' in a list
#' @return if `keep = FALSE`, a data.frame summarizing the search with query name, maching
#' subject, percentage of their identity, alignment length, and number of mismatches.
#' Where subject wasn't identified, values are `NA`.
#' If `keep = TRUE`, a list containing the above data.frame as well as results from the forward
#' and backward searches, together with the best hits in both directions.
#'
#' @export
rblast = function(
    x,
    query,
    subject,
    type = c("blastn", "blastp"),
    query_db = NULL,
    subject_db = NULL,
    keep = FALSE
    ){

    type = match.arg(type)

    x = strfsplit(x, " ", fixed = TRUE)[[1]]

    query_names_match = match_names(x, names(query))

    db_type = switch(type,
        "blastn" = "nucl",
        "blastp" = "prot"
        )

    if(is.null(query_db)){
        query_db = file.path(tempdir(), "query")
        make_blast_db(query, out = query_db, type = db_type)
        }

    if(is.null(subject_db)){
        subject_db = file.path(tempdir(), "subject")
        make_blast_db(subject, out = subject_db, type = db_type)
        }

    # first blast
    forward = blast(
        query[query_names_match],
        db = subject_db,
        outfmt = 6,
        type = type,
        args = "-max_target_seqs 5"
        )
    forward_best = forward |> split(~ query) |> lapply(utils::head, 1) |> do.call(what = rbind)

    subject_names = names(subject)
    subject_names_match = match_names(forward_best$subject, subject_names)

    # second blast
    backward = blast(
        subject[subject_names_match],
        db = query_db, outfmt = 6,
        type = type,
        args = "-max_target_seqs 5"
        )
    backward_best = backward |> split(~ query) |> lapply(utils::head, 1) |> do.call(what = rbind)
    backward_best = backward_best[forward_best$subject, ]

    res = data.frame(
        "query" = forward_best$query,
        "subject" = forward_best$subject,
        "perc_identity" = forward_best$perc_identity,
        "alignment_length" = forward_best$alignment_length,
        "mismatches" = forward_best$mismatches,
        "annotation" = subject_names[subject_names_match]
        )

    # remove what is actually not a match
    res[forward_best$query != backward_best$subject, -1] = NA
    res = res[match(x, res$query),]
    rownames(res) = NULL
    res$query = x

    if(keep){
        res = list(
            "rblast" = res,
            "forward" = forward,
            "forward_best" = forward_best,
            "backward" = backward,
            "backward_best" = backward_best
            )
        }

    res
    }


match_names = function(x, where){
    where = strfsplit(where, " ", fixed = TRUE)[[1]]
    y = lapply(x, \(y) which(y == where)) |> stats::setNames(x)

    # sanity check
    if(any(lengths(y) != 1))
        stop("Non-unique match of sequence names! This shouldn't happen.")

    unlist(y)
    }
