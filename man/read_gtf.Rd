% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtf.r
\name{read_gtf}
\alias{read_gtf}
\alias{parse_attributes}
\title{Read gtf file}
\usage{
read_gtf(x, feature = NULL, attributes = FALSE)

parse_attributes(x)
}
\arguments{
\item{x}{For \code{read_gtf()} A file in the gtf2 format, for \code{parse_attributes()} the attribute
column see details}

\item{feature}{Filter the gtf file according to a feature of interest, such as gene, CDS,
or exon. For instance, if you are interested in proteins, they are typically carried by the
CDS feature.}

\item{attributes}{Parse attributes using the `parse_attributes() function, note that this
might take a substantial time.}
}
\value{
For \code{read_gtf()} a data.frame with a parsed content of the gtf file with
the number of rows equal to the number of annotated lines in the gtf file, excluding comments.
If \code{attributes = FALSE}, the \code{data.frame} will contain exactly 9 columns (see details).
If \code{attributes = TRUE}, the number of columns varies depending on the annotation,
see details.
For \code{parse_attributes()}, a data.frame with parsed attributes with the number of rows
equivalent to the length of the input vector.
}
\description{
Read the gtf file and optionally parse attributes.
}
\details{
The GTF2.2 file format is an extension of the older GFF file format with additional
information. Typically, GTF2.2 has 9 tab-separated columns: chromosome, source, feature,
start, end, score, strand, frame, and the attribute column. The attribute column is optional,
but in practice almost always present (at least the author of this comment have not seen
a gtf file without it).

If the \code{attributes = TRUE}, the attribute field is parsed and the attribute column removed.
Currently, parsing attributes is performed fully in R and is rather slow, filtering by
feature and caching is recommended. The field consist of semicolon-separated fields,
their number can vary depending on the available annotation information.
The field themselves consist of two space-separated elements, first one defines
the name of the field with the second being a double-quoted value, and there can be
multiple fields with the same name.

The \code{parse_attributes()} function parses these fields, creating a \code{data.frame} with
all shared fields, if multiple fields share the same name, they are merged using semicolon.
}
\seealso{
Definition of GTF 2.2: http://mblab.wustl.edu/GTF22.html
}
