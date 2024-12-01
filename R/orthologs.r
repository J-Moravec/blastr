#' Find orthologs based on RBH protein search.
#'
#' This is a user-friendly function that downloads query and subject protein sequences,
#' and then performs a BLAST reciprocal best hit (RBH) search to identify ortholog sequences.
#' The input query and output subject genes are matched to proteins using the downloaded gtf
#' annotation files.
#'
#' This is an example function that employs many parts of the `blastr` package.
#' First, the query and subject genomes are downloaded from the NCBI using the [ncbi()] function.
#' The input gene names are matched to protein names using the query `gtf` file.
#' Then, RBH is run to identify matched subject protein sequences, which are then translated
#' to gene names using the subject `gtf` file.
#'
#' Since parsing the `gtf` file might take some time, for repeated runs the parsed
#' `gtf` file can be cached.
#'
#' @param genes Gene names for which ortholog sequences should be found.
#' @param query,subject NCBI genome ID for query and subject respectively.
#' @param dir A directory where the NCBI genome data are downloaded,
#' the default is the current working directory.
#' @param cache A path to cache file. The default is ".cache/annotation.rds".
#' Set to `cache = NULL` to disable caching.
#' @param nthreads The number of threads parameter passed to NCBI BLAST+.
#' The default is `1`.
#' @param keep Keep the forward and backward (reciprocal) blast searches,
#' @param type Run the `RBH` or `RBSH` method, default is `RBH`
#' these are returned together with normal output in a named list
#' @param ... Other arguments passed to `rblast` or `rsblast` respectively.
#'
#' @return A data.frame of with `length(genes)` rows and 5 columns:
#' `query_gene`, `query_protein`, `subject_gene`, `subject_protein`, and `subject_annotation`.
#' Where query gene couldn't be matched to query protein, or the query protein didn't rblast
#' any subject protein, the values are `NA_character_`.
#'
#' @export
#' @seealso
#' [blastr::ncbi()] for downloading the genomes from the NCBI genome database,
#' [blastr::read_fasta()] and [blastr::read_gtf()] for reading the `fasta` and `gtf` files
#' respectively,
#' [blastr::filecache()] for a simple on-disk caching,
#' [blastr::rblast()] and [blastr::rsblast()] for the implementation of RBH and RBSH
orthologs = function(
    genes, query, subject,
    dir = ".", cache = ".cache/annotation.rds",
    nthreads = 1, keep = FALSE, type = c("RBH", "RBSH"),
    ...
    ){
    faa = \(x) file.path(dir, paste0(x, ".faa.gz"))
    gtf = \(x) file.path(dir, paste0(x, ".gtf.gz"))

    type = match.arg(type)
    blastr::ncbi(query, formats = c("PROT_FASTA", "GENOME_GTF"), dir = dir)
    blastr::ncbi(subject, formats = c("PROT_FASTA", "GENOME_GTF"), dir = dir)

    query_seq = faa(query) |> blastr::read_fasta()
    subject_seq = faa(subject) |> blastr::read_fasta()

    if(is.null(cache)){
        query_annotation = read_gtf(
            gtf(query),
            feature = "CDS",
            attributes = TRUE
            )[c("gene_id", "protein_id")]
        subject_annotation = read_gtf(
            gtf(subject),
            feature = "CDS",
            attributes = TRUE
            )[c("gene_id", "protein_id")]
        } else {
        query_annotation = filecache(
            cache,
            read_gtf, gtf(query), feature = "CDS", attributes = TRUE
            )[c("gene_id", "protein_id")]
        subject_annotation = filecache(
            cache,
            read_gtf, gtf(subject), feature = "CDS", attributes = TRUE
            )[c("gene_id", "protein_id")]
        }

    map = data.frame(
        "query_protein" = strfsplit(names(query_seq), " ", fixed = TRUE)[[1]]
        )
    map[["query_gene"]] = query_annotation[
        match(map[["query_protein"]], query_annotation[["protein_id"]]),
        "gene_id"
        ]
    map = merge(data.frame("query_gene" = genes), map, all.x = TRUE)

    res = switch(type,
        "RBH" = rblast(
            map$query_protein |> na.rm(),
            query_seq,
            subject_seq,
            type = "blastp",
            args = paste0("-num_threads ", nthreads),
            keep = keep,
            ...
            ),
        "RBSH" = rsblast(
            map$query_protein |> na.rm(),
            query_seq,
            subject_seq,
            type = "blastp",
            args = paste0("-num_threads ", nthreads),
            keep = keep,
            ...
            )
        )

    if(keep){
        details = res[-1]
        res = res[[1]]
        }

    res = res[c("query", "subject", "annotation", "perc_identity", "alignment_length")]
    names(res) = c("query_protein", "subject_protein", "subject_annotation", "subject_identity", "subject_alength")

    res[["subject_gene"]] = subject_annotation[
        match(res[["subject_protein"]], subject_annotation[["protein_id"]]),
        "gene_id"
        ]
    res[["subject_annotation"]] = res[["subject_annotation"]] |>
        strfsplit(" ", fixed = TRUE) |>
        getElement(2) |>
        sub(pattern = " \\[[A-Za-z ]*\\]$", replacement = "")

    map = merge(map, res, all = TRUE)[c(
        "query_gene", "query_protein",
        "subject_gene", "subject_protein", "subject_annotation", "subject_identity", "subject_alength"
        )]

    if(keep){
        map = list("search" = map, "details" = details)
        }

    map
    }


na.rm = function(x){
    x[!is.na(x)]
    }
