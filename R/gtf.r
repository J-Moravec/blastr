#' Read gtf file
#'
#' Read the gtf file and optionally parse attributes.
#'
#' The GTF2.2 file format is an extension of the older GFF file format with additional
#' information. Typically, GTF2.2 has 9 tab-separated columns: chromosome, source, feature,
#' start, end, score, strand, frame, and the attribute column. The attribute column is optional,
#' but in practice almost always present (at least the author of this comment have not seen
#' a gtf file without it).
#'
#' If the `attributes = TRUE`, the attribute field is parsed and the attribute column removed.
#' Currently, parsing attributes is performed fully in R and is rather slow, filtering by
#' feature and caching is recommended. The field consist of semicolon-separated fields,
#' their number can vary depending on the available annotation information.
#' The field themselves consist of two space-separated elements, first one defines
#' the name of the field with the second being a double-quoted value, and there can be
#' multiple fields with the same name.
#'
#' The `parse_attributes()` function parses these fields, creating a `data.frame` with
#' all shared fields, if multiple fields share the same name, they are merged using semicolon.
#'
#' @param x For `read_gtf()` A file in the gtf2 format, for `parse_attributes()` the attribute
#' column see details
#' @param feature Filter the gtf file according to a feature of interest, such as gene, CDS,
#' or exon. For instance, if you are interested in proteins, they are typically carried by the
#' CDS feature.
#' @param attributes Parse attributes using the `parse_attributes() function, note that this
#' might take a substantial time.
#' @return For `read_gtf()` a data.frame with a parsed content of the gtf file with
#' the number of rows equal to the number of annotated lines in the gtf file, excluding comments.
#' If `attributes = FALSE`, the `data.frame` will contain exactly 9 columns (see details).
#' If `attributes = TRUE`, the number of columns varies depending on the annotation,
#' see details.
#' For `parse_attributes()`, a data.frame with parsed attributes with the number of rows
#' equivalent to the length of the input vector.
#'
#' @export
#'
#' @seealso
#' Definition of GTF 2.2: http://mblab.wustl.edu/GTF22.html
read_gtf = function(x, feature = NULL, attributes = FALSE){
    # gtf contains comments within fields as well
    y = readLines(x)
    y = y |> grep(pattern = "^#", value = TRUE, invert = TRUE)
    y = scan(text = y, sep = "\t", quote = "", what = rep(list(character()), 9), quiet = TRUE)
    names(y) = c("chr", "source", "feature", "start", "end",
                 "score", "strand", "frame", "attribute")
    y = data.frame(y)

    if(!is.null(feature))
        y = y[y$feature %in% feature,]

    if(attributes)
        y = cbind(y[-9], y$attribute |> parse_attributes())

    y
    }


#' @rdname read_gtf
#' @export
parse_attributes = function(x){
    y = strsplit(x, "\"; ", fixed = TRUE)
    y = lapply(y, parse_attribute_fields)
    y = rbindfill(y)
    y[] = lapply(y, gsub, pattern = "\"", replacement = "")
    y
    }


parse_attribute_fields = function(x, sep = " "){
    s = strfsplit(x, pattern = sep, fixed = TRUE)
    aggregate(s[[2]], s[[1]])
    }


rbindfill = function(x){
    fill = function(x, names){
        x[setdiff(names, names(x))] = NA
        x
        }

    nm = lapply(x, names) |> unlist() |> unique()
    x = lapply(x, fill, names = nm)
    do.call(rbind.data.frame, x)
    }


aggregate = function(x, by){
    if(length(x) != length(by))
        stop("x and by must have the same length!")

    u = unique(by)
    lapply(u, \(y) paste0(x[by == y], collapse = "; ")) |> stats::setNames(u)
    }
