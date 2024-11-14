# Ported and modified from:
# https://github.com/J-Moravec/ncbi.r

basename_sans_ext = function(x){
    tools::file_path_sans_ext(basename(x))
    }

"%nin%" = Negate("%in%")

ncbi_api_url =  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"


#' NCBI genomic formats
#'
#' Functions to work with NCBI genomic file formats.
#'
#' The NCBI datasets API define file formats that are available to download together
#' with the genome. These formats carry additional information, such as coding sequences,
#' annotation, or translated protein sequences. However, not all formats are available for
#' all genomes. In the API documentation, they are refered to as annotation type, here
#' we refer to them as formats, but since most of these files are in the `fasta` format
#' and only `.gtf`, `.gff`, and `.gbff` files are considered annotation, neither term
#' is accurate.
#'
#' The `ncbi_formats` object is a named vector of all possible formats.
#' The names are used in the NCBI JSON output, while the values are used in the API input,
#' as well as in the other ncbi functions.
#'
#' The `ncbi_available_formats()` makes a NCBI API request to find what formats
#' are actually available for the requested genome.
#'
#' The `ncbi_format_mapping_table()` is then used to translate the file names between
#' unflattened and flattened archive.
#'
#' @param id NCBI genome accession id
#' @param assembly_name an assembly name obtained from parsing `assembly_data_report.jsonl`.
#' Normally, the assembly has the same name as the genome accession id, but in some cases it
#' might differ.
#' @param formats **optional** requested formats to display
#' @return In case of `ncbi_available_formats()`, a character vector with a subset of values
#' from `ncbi_formats`. The `ncbi_format_mapping_table()` returns a data.frame with original
#' and new file names.
#'
#' @export
ncbi_formats = c(
    "all_genomic_fasta" = "GENOME_FASTA",
    "genome_gff" = "GENOME_GFF",
    "genome_gtf" = "GENOME_GTF",
    "genome_gbff" = "GENOME_GBFF",
    "rna_fasta" = "RNA_FASTA",
    "cds_fasta" = "CDS_FASTA",
    "prot_fasta" = "PROT_FASTA",
    "sequence_report" = "SEQUENCE_REPORT"
    )


#' @rdname ncbi_formats
#' @export
ncbi_available_formats = function(id){
    url = paste(ncbi_api_url, id, "download_summary", sep = "/") |> url()
    genome_summary = readLines(url, warn = FALSE) # missing final line
    # transforming complex json into a simpler structure by ... deleting stuff
    genome_summary = genome_summary |>
        strsplit(split = "[{},]") |>
        getElement(1) |> trimws() |>
        gsub(pattern = "\"", replacement = "") |>
        gsub(pattern = ":.*", replacement = "")
    # now "file_count" is just under file name
    available_formats = genome_summary[which(genome_summary == "file_count") - 1]
    # now we need to transform the file names into formats
    ncbi_formats[available_formats]
    }


#' @rdname ncbi_formats
#' @export
ncbi_format_mapping_table = function(id, assembly_name, formats = NULL){
    genome_fasta = paste(id, assembly_name, "genomic.fna", sep = "_")

    id_ext = function(ext){
        paste0(id, ".", ext)
        }

    mapping = list(
        c("GENOME_FASTA", genome_fasta, id_ext("fna")),
        c("GENOME_GFF", "genomic.gff", id_ext("gff")),
        c("GENOME_GTF", "genomic.gtf", id_ext("gtf")),
        c("RNA_FASTA", "rna.fna", id_ext("rna.fna")),
        c("CDS_FASTA", "cds_from_genomic.fna", id_ext("cds.fna")),
        c("PROT_FASTA", "protein.faa", id_ext("faa")),
        c("SEQUENCE_REPORT", "sequence_report.jsonl", id_ext("jsonl"))
        ) |> do.call(what = rbind.data.frame) |> stats::setNames(c("formats", "old", "new"))
    mapping[["old"]] = file.path("ncbi_dataset", "data", id, mapping[["old"]])

    if(!is.null(formats))
        mapping = mapping[match(formats, mapping[["formats"]], nomatch = 0),]

    mapping
    }


#' Download, extract, and flatten the NCBI genome archives
#'
#' Functions to download, extract, and flatten the NCBI genome archives.
#'
#' The NCBI genome archive allows you to download reference genomes through a web interface
#' or through REST API in a nice zip archive. There is only one issue, it comes in a deeply
#' nested structure with a lot of ancillary files.
#'
#' The purpose of these functions is to download the NCBI archive, unpack it, extract
#' the files of interest, and rename them with a unified format with `{accession}.{extension}`.
#'
#' The `ncbi()` downloads, extract, and flattens the archive. If `formats = NULL`,
#' it will perform an additional NCBI API request to determine which formats are available
#' for the requested genome.
#'
#' `ncbi_download()` downloads a zip archive from the NCBI genome website.
#'
#' `ncbi_extract()` is a simple wrapper around [unzip()] that extracts the NCBI genome archive.
#'
#' `ncbi_flatten()` uses `ncbi_extract()` internally to flatten the archive. Flatteining
#' in here means extracting the requested files from a deeply nested structure and returning
#' a simple list of consistently named files.
#'
#' @param id NCBI accession id of requested genome
#' @param formats requested file formats, such as `GENOME_FASTA` for fasta sequences
#' or `GENOME_GTF` for the annotation in the `.gtf` format.
#' For the `ncbi()` function, if `formats = NULL`, the available formats for the accession `id`
#' are obtained and all of them are downloaded. This might duplicate some information,
#' for instance the Genbank Flat File Format `gbff` (`GENOME_GBFF`) includes sequences
#' and annotation information, if you ar ealready downloading `GENOME_FASTA` and `GENOME_GTF`,
#' you might consider to leave out this format to save space and bandwidth.
#' `ncbi_flatten()` need to have formats specified and defaults to `c("GENOME_FASTA", "GENOME_GTF")`.
#' @param dir **optional** directory where the files will be downloaded. For `ncbi()`,
#' `ncbi_flatten()`, and `ncbi_extract()`, the default folder is a folder with the same
#' name as the accession `id`. The `ncbi_download()` defaults to the current working directory.
#' @param compress **optional** compress the extracted files using [gzip()].
#' @param overwrite **optional** overwrite existing files.
#' @param keep **optional** keep the downloaded archive
ncbi = function(
    id,
    formats = NULL,
    dir = NULL,
    compress = TRUE,
    overwrite = FALSE,
    keep = FALSE
    ){
    if(is.null(formats))
        formats = ncbi_available_formats(id)

    formats = match.arg(formats, ncbi_formats, several.ok = TRUE)

    if(is.null(dir))
        dir = id

    targets = file.path(dir, ncbi_format_mapping_table(id, "phony", formats)$new)
    if(compress)
        targets = paste0(targets, ".gz")

    redo = overwrite || !all(file.exists(targets))

    if(!redo)
        return(invisible(targets))

    var = id
    var = ncbi_download(var, dir = dir, formats = formats, overwrite = overwrite)
    var = ncbi_flatten(var, dir = dir, formats = formats, keep = keep, overwrite = overwrite)
    if(compress)
        var = sapply(var, gzip, overwrite = overwrite)

    invisible(var)
    }



#' @rdname ncbi
#' @param quiet **optional** print diagnostic information during download,
#' passed to [download.file()]
#' @param timeout **optional** maximum allowed time for the connection in seconds,
#' passed to [download.file()]
#' @export
ncbi_download = function(
    id,
    dir = ".",
    formats = c("GENOME_FASTA", "GENOME_GTF"),
    overwrite = FALSE,
    quiet = TRUE,
    timeout = 3600
    ){
    dest = file.path(dir, paste0(id, ".zip"))

    if(file.exists(dest) && !overwrite)
        return(invisible(dest))

    formats = match.arg(formats, ncbi_formats, several.ok = TRUE)

    formats = paste0(formats, collapse = ",")

    url = paste0(
        ncbi_api_url,
        id, "/download?include_annotation_type=",
        formats
        )

    # Remove file if something stopped the download
    ok = FALSE
    on.exit(
        if(!ok){
            file.remove(dest)
            stop("Download of '", id, "' from NCBI failed.", call. = FALSE)
            },
        add = TRUE
        )

    op = options("timeout" = timeout)
    err = utils::download.file(url, dest, quiet = quiet)
    ok = TRUE
    options(op)

    if(err != 0){
        file.remove(dest)
        stop(paste0("Download of '", id, "' from NCBI failed."))
        }

    # If invalid ID is provided, the download progresses normally but the file is malformed.
    # The usual size of such file seems to be 855 bytes.
    if(file.size(dest) < 1000){
        file.remove(dest)
        stop(paste0("The requested ID '", id, "' doesn't exist."))
        }

    invisible(dest)
    }


#' @rdname ncbi
#' @param archive a path to the downloaded NCBI genomic archive
#' @export
ncbi_flatten = function(
    archive, dir = NULL,
    formats = c("GENOME_FASTA", "GENOME_GTF"),
    keep = FALSE, overwrite = FALSE
    ){
    formats = match.arg(formats, ncbi_formats, several.ok = TRUE)
    prefix = basename_sans_ext(archive)

    if(is.null(dir))
        dir = prefix

    # early exit
    new_files = file.path(dir, ncbi_format_mapping_table(prefix, "phony", formats)$new)
    if(all(file.exists(new_files)) && !overwrite)
        return(invisible(new_files))

    y = ncbi_extract(archive, dir = dir, keep = keep)
    assembly_name = get_assembly_name(
        file.path(y, "ncbi_dataset", "data", "assembly_data_report.jsonl")
        )
    mapping = ncbi_format_mapping_table(prefix, assembly_name, formats)

    # safety check
    files = list.files(y, recursive = TRUE)
    if(!all(mapping$old %in% files)){
        missing = mapping$formats[mapping$old %nin% files]
        stop("Some of the requested formats are not in the NCBI archive:\n", toString(missing))
        }

    # move and remove
    old_files = file.path(y, mapping$old)
    file.rename(old_files, new_files) # no copying involved, should be faster

    # always remove the folder
    if(TRUE)
        unlink(y, recursive = TRUE, force = TRUE)

    invisible(new_files)
    }


#' @rdname ncbi
#' @export
ncbi_extract = function(archive, dir = NULL, keep = FALSE){
    ext = tools::file_ext(archive)
    prefix = basename_sans_ext(archive)

    if(is.null(dir))
        dir = prefix

    if(ext != "zip")
        stop(paste0("Unrecognised extension '", ext, "'. Only zip archives are supported.")) 

    utils::unzip(archive, exdir = dir, unzip = "unzip")

    if(!keep)
        file.remove(archive)

    invisible(dir)
    }


get_assembly_name = function(x){
    y = readLines(x)
    regmatches(y, regexpr("assemblyName\":\"([^\"]*)", y)) |> substring(16)
    }
