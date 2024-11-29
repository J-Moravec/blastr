#' Reciprocal BLAST
#'
#' Perform BLAST reciprocal best hit (RBH) search.
#'
#' BLAST reciprocal best hit search will first perform a forward search, take the best hits,
#' and then performs a backward search using these best hits. If the sequence used for
#' the forward search are identified as a best hit in the backward search, the matched
#' subject are reported.
#'
#' The function `rblast` implements the standard reciprocal blast.
#' The function `rsblast` implements an extended version of reciprocal blast called
#' Reciprocal Best Similar Hit blast. In this method, the forward blast is performed as normal,
#' but in the backward (reciprocal) blast instead of taking only the best hit,
#' all hits that are similar enough (in terms of `bit_score`) to the best hit are considered.
#' This should outperform RBH in a situation where many orthologs exist.
#'
#' Internally `rblast` is implemented within `rsblast` when `p = NULL`.
#'
#' @param x Names or ids of the sequences to be blasted, id is the first space-separated part
#' of a sequence name and is intepreted by the NCBI BLAST+ as a query name.
#' @param query,subject A character vector of sequences of the query and target organisms,
#' such as one from [read_fasta()].
#' @param type a type of blast, either nucleotide `blastn` or protein `blastp` blast
#' @param query_db,subject_db pre-existing query and subject blast dabatases.
#' If not provided, they are created in a R's tempdir from the sequences, see [make_blast_db()].
#' @param args a character vector of additional arguments passed to blast.
#' @param keep Keep the forward and backward (reciprocal) blast searches,
#' these are returned together with normal output in a named list
#' @param sort_by A character vector of column names to sort by when chosing the best hit.
#' By default, the values are sorted by `evalue`. Other options to consider might be `bit_score`, or `perc_identity`,
#' but `evalue` should already provide a comprehensive metric for the sequence match.
#' @param decreasing A logical vector ideally of the same length as `sort_by`.
#' Whether the columns specified by `sort_by` should be sorted in decreasing order.
#' Default for `evalue` is `FALSE` since smaller `evalue` is more significant.
#' @param n The number of subject returned for each blast search (`-max_target_seqs`), defaults to 5.
#' @param p The proportion of `bit_score` of the best hit within which other hits will be considered.
#' If `p = NULL`, standard RBH is run instead.
#' @return if `keep = FALSE`, a data.frame summarizing the search with query name, maching
#' subject, percentage of their identity, alignment length, and number of mismatches.
#' Where subject wasn't identified, values are `NA`.
#' If `keep = TRUE`, a list containing the above data.frame as well as results from the forward
#' and backward searches, together with the best hits in both directions.
#'
#' @export
rblast = function(
    x, query, subject,
    type = c("blastn", "blastp"),
    query_db = NULL, subject_db = NULL,
    args = NULL, keep = FALSE,
    sort_by = "evalue", decreasing = FALSE,
    n = 5
    ){
    rsblast(x, query, subject, type, query_db, subject_db, args, keep, sort_by, decreasing, p = NULL, n)
    }


#' @rdname rblast
#' @export
rsblast = function(
    x,
    query,
    subject,
    type = c("blastn", "blastp"),
    query_db = NULL,
    subject_db = NULL,
    args = NULL,
    keep = FALSE,
    sort_by = "evalue",
    decreasing = FALSE,
    p = 0.1,
    n = 5
    ){
    type = match.arg(type)
    args = c(
        args,
        paste0("-max_target_seqs ", as.integer(n[1]))
        )

    x = strfsplit(x, " ", fixed = TRUE)[[1]]

    query_names_match = match_names(x, names(query))

    db_type = switch(type,
        "blastn" = "nucl",
        "blastp" = "prot"
        )

    if(is.null(query_db))
        query_db = make_blast_db_tmp(query, "query.db", db_type)
    if(is.null(subject_db))
        subject_db = make_blast_db_tmp(subject, "subject.db", db_type) 

    forward = blast(
        query[query_names_match],
        db = subject_db,
        outfmt = 6,
        type = type,
        args = args
        )
    forward_best = forward |>
        split(~ query) |>
        lapply(best_sorted_hit, sort_by = sort_by, decreasing = decreasing) |>
        do.call(what = rbind)

    subject_names = names(subject)
    subject_names_match = match_names(forward_best$subject, subject_names)

    backward = blast(
        subject[subject_names_match],
        db = query_db, outfmt = 6,
        type = type,
        args = args
        )

    if(is.null(p)){
        backward_best = backward |>
            split(~ query) |>
            lapply(best_sorted_hit, sort_by = sort_by, decreasing = decreasing) |>
            do.call(what = rbind)
        backward_best = backward_best[forward_best$subject, ]
        is_ortholog = forward_best$query == backward_best$subject
        } else {
        backward_best = backward |>
            split(~ query) |>
            lapply(best_sorted_bit_scores, sort_by = sort_by, decreasing = decreasing, p = p)
        backward_best = backward_best[forward_best$subject]
        backward_best_subject = lapply(backward_best, getElement, "subject")
        backward_best = do.call(what = rbind, backward_best)
        is_ortholog = mapply(`%in%`, forward_best$query, backward_best_subject, SIMPLIFY = TRUE)
        }

    res = data.frame(
        "query" = forward_best$query,
        "subject" = forward_best$subject,
        "perc_identity" = forward_best$perc_identity,
        "alignment_length" = forward_best$alignment_length,
        "mismatches" = forward_best$mismatches,
        "annotation" = subject_names[subject_names_match]
        )

    res[!is_ortholog, -1] = NA
    res = res[match(x, res$query),]
    rownames(res) = NULL
    res$query = x

    if(keep){
        res = list(
            "search" = res,
            "forward" = forward,
            "forward_best" = forward_best,
            "backward" = backward,
            "backward_best" = backward_best
            )
        }

    res
    }


best_sorted_hit = function(x, sort_by, decreasing){
    sort_by(x, x[sort_by], decreasing = decreasing) |> utils::head(1)
    }


best_sorted_bit_scores = function(x, sort_by, decreasing, p){
    x = sort_by(x, x[sort_by], decreasing = decreasing)
    threshold = max(x[["bit_score"]], na.rm = TRUE) * (1 - p)
    x[which(x[["bit_score"]] > threshold),]
    }


make_blast_db_tmp = function(x, name = "db", type = c("nucl", "prot")){
    type = match.arg(type)

    file = file.path(tempdir(), name)
    make_blast_db(x, out = file, type = type)

    file
    }
