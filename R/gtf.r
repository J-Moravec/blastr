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
    y = scan(x, sep = "\t", quote = "", what = rep(list(character()), 9),
             comment.char = "", fill = TRUE, quiet = TRUE)
    names(y) = c("chr", "source", "feature", "start", "end",
                 "score", "strand", "frame", "attribute")
    y = data.frame(y)
    y = y[!startsWith(y[[1]], "#"),]

    if(!is.null(feature))
        y = y[y$feature %in% feature,]

    if(attributes)
        y = cbind(y[-9], y$attribute |> parse_attributes())

    rownames(y) = NULL
    y
    }


#' @rdname read_gtf
#' @export
parse_attributes = function(x){
    # trailing "; " caused incorrect parsing
    y = sub(x, pattern = ";[ ]*$", replacement = "")
    y = strsplit(y, "\"; ", fixed = TRUE)
    for(i in seq_along(y)){
        z = y[[i]] |> strfsplit(pattern = " ", fixed = TRUE)
        y[[i]] = aggregate(z[[2]], z[[1]])
        }
    y = rbindfill(y)
    y = gsub(y, pattern = "\"", replacement = "", fixed = TRUE)
    as.data.frame(y)
    }


rbindfill = function(x){
   fill = function(x, names){
        y = rep(NA_character_, length(names))
        names(y) = names
        y[names(x)] = x
        y
        }

    nm = lapply(x, names) |> unlist() |> unique()
    vapply(x, fill, character(length(nm)), names = nm) |> t()
    }


aggregate = function(x, by){
    u = unique(by)
    y = vapply(
        u,
        \(y){
            z = x[by == y]
            if(length(z) == 1) z else paste0(z, collapse = "; ")
            },
        character(1)
        )
    names(y) = u
    y
    }
