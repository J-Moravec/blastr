#' Make BLAST database
#'
#' Make BLAST database out of fasta file.
#'
#' @param x fasta file to make a database from, will be unzipped if required
#' @param out an output path and name of the database, such as `path/name.db`
#' @param type the type of database
#' @param title **optional** BLAST database title, defaults to the input file name
#' @param verbose **optional** prints diagnostic information
#'
#' @export
make_blast_db = function(x, out, type = c("nucl", "prot"), title = NULL, verbose = FALSE){
    if(is_sequences(x))
        x = write_fasta_temp(x)

    if(identical(tools::file_ext(x), "gz")){
        y = tools::file_path_sans_ext(x) |> basename()
        x = gunzip(x, file.path(tempdir(), y), keep = TRUE)
        }

    command = "makeblastdb"
    check_prog(command)

    type = match.arg(type)

    args = c(
        "-in", x,
        "-out", out,
        "-dbtype", type,
        if(!is.null(title)) c("-title", title)
        )

    stdout = TRUE
    if(verbose)
        stdout = ""

    system2(command, args, stdout = stdout) |> invisible()
    }
