#' Compress or decompress files
#'
#' Compress or decompress input file and write the encoded byte-stream into output file.
#'
#' The compression or decompression depends entirely on the passed `input_con` and `output_con`
#' functions. If `input_con` is `file()` and `output_con` a `gzfile()`, then the function performs
#' standard gzip compression.
#'
#' @param input,output a file for input/output
#' @param input_con,output_con a function to create input/output connections
#' @param overwrite overwrite output if it exists
#' @param keep preserve input
#'
#' @seealso
#' Connections like [base::file()], [base::gzfile()], [base::bzfile()], and [base::xzfile()]
#' that can be used by this function to work with gz, bz, or xz file compressions.
#' [gzip()] and [gunzip()] which are simple user-friendly wrappers around this function.
#'
#' @export
compress_or_decompress = function(
    input, output,
    input_con, output_con,
    overwrite = FALSE, keep = FALSE
    ){
    if(input == output)
        stop("Input and output paths are identical.")

    if(file.exists(output) && !overwrite)
        return(invisible(output))

    if(file.exists(output))
        file.remove(output)

    completed = FALSE
    input_connection = input_con(input, open = "rb")
    on.exit(expr = {
        close(input_connection)
        if(completed && !keep) file.remove(input)
        }, add = TRUE)

    output_connection = output_con(output, open = "wb")
    on.exit(expr = {
        close(output_connection)
        if(!completed) file.remove(output)
        }, add = TRUE)

    repeat {
        b = readBin(input_connection, what = raw(0L), size = 1L, n = 1e7)
        if (length(b) == 0L) break
        writeBin(b, con = output_connection, size = 1L)
        }

    completed = TRUE
    invisible(output)
    }


#' Compress and decompress using gzip
#'
#' This is a native implementation of gzip using the `base::gzfile()` connection.
#'
#' R can natively read and write in the `gz` format, but miss functions compress
#' and decompress files on disk. In this implementation, files are not being
#' read instead the compression happens on a chunked bytestream, so even
#' large files can be compressed or decompressed without exhausting RAM.
#'
#' The interface tries to mimic the unix tools `gzip` and `gunzip`.
#' If `out` is not provided, the output file is derived from the input file, either
#' appending (for `gzip()`) or removing (for `gunzip()`) the `.gz` extension.
#' Once the file is successfully compressed/decompressed, the input file is removed
#' unless `keep = TRUE` is specified.
#' If for some reason the process fails, broken files are removed.
#'
#' @param x file to be compressed or decompressed using the `gz` algorithm
#' @param out **optional** output file, see details
#' @param overwrite **optional** overwrite existing file
#' @param keep **optional** keep the original file
#' @return the output path is silently returned
#'
#' @examples
#' # Compression and decompression
#' test = tempfile()
#' writeLines("This is a test!", test)
#' test.gz = gzip(test)
#' file.exists(test) # no longer exists!
#' test = gunzip(test.gz)
#' file.exists(test.gz) # no longer exists!
#' readLines(test) # "This is a test!"
#'
#' @seealso
#' [base::gzfile()] for the gz connection
#' [compress_or_decompress()] for the wokrhorse of `gzip` and `gunzip`
#'
#' @export
gzip = function(x, out = NULL, overwrite = FALSE, keep = FALSE){
    if(is.null(out))
        out = paste0(x, ".gz")

    compress_or_decompress(
        x, out,
        input_con = file,
        output_con = gzfile,
        overwrite = overwrite,
        keep = keep
        )
    }


#' @rdname gzip
#' @export
gunzip = function(x, out = NULL, overwrite = FALSE, keep = FALSE){
    if(is.null(out))
        out = tools::file_path_sans_ext(x)

    compress_or_decompress(
        x, out,
        input_con = gzfile,
        output_con = file,
        overwrite = overwrite,
        keep = keep
        )
    }
