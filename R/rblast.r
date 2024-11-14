#' Reciprocal BLAST
#'
#' Perform BLAST reciprocal best hit (RBH) search.
#'
#' @param x names of sequences to be blasted
#' @param query,subject sequences of the query and target organisms
#' @param query_db,subject_db **optional** pre-existing database made from query and subject
#' sequences respectively
#' @param type a type of blast, either nucleotide `blastn` or protein `blastp` blast
#' @param query_name,subject_name **optional** names for the query and subject databases
#' @param dir = "." **optional** wokring directory
#'
#rblast = function(){
#    x,
#    query,
#    subject,
#    query_db = NULL,
#    subject_db = NULL,
#    type = c("blastn", "blastp")
#    query_name = NULL,
#    subject_name = NULL,
#    dir = "."
#    }

#    type = match.arg(type)

#    if(is.null(query_name))
#        query_name = "query"

#    if(is.null(subject_name))
#        subject_name = "subject"

#    if(is.null(query_db)){
#        make_blast_db(query)
#        }

#    if(is.null(subject_db)){
#        make_blast_db(subject)
#        }

#    # first blast
#    blast

#    # second blast
#    blast
