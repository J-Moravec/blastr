# blastr: mapping to NCBI BLAST

`blastr` is a lightweight no-dependency package providing a mapping to the [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), a local version of [BLAST](https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi).

Outside of `blast()` and its wrappers `blastn()`, `blastp()`, etc. for blasting,
`make_blast_db()` for creating BLAST database, and `blast_formatter()` for formatting BLAST output,
the package `blastr` also provide functions for downloading reference genomes from NCBI
(taken from [ncbi.r](https://github.com/J-Moravec/ncbi.r), and a functions to perform a large
number of blast searches at the same time, to improve workflows such as identification of
ortholog genes.

## Installation

Use your favorite of `pak`, `remotes`, or `base::install.packages()`

## Usage
