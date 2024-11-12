#' Run BLAST query
#'
#' Perform BLAST search
#'
#' @param query a fasta file or a fasta object
#' @param subject **optional** fasta file or a fasta object against which to blast, see details
#' @param db **optional** a database created from [make_blast_db()], see details
#' @param type **optional** type of BLAST to run, one of `blastn`, `blastp`, `blastx`,
#' `tblastn` or `tblastx`.
#' @param out output path
#' @param outfmt output format
#' @return dunno
blast = function(query, subject = NULL, db = NULL, type = "blastn", out, outfmt){
    #TODO
    stop("Unimplemented")

    }
