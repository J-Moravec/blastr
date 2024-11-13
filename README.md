# blastr: mapping to NCBI BLAST

`blastr` is a lightweight no-dependency package providing a binding to the [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), a local version of [BLAST](https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi).

Outside of `blast()` and its wrappers `blastn()`, `blastp()`, etc. for blasting,
`make_blast_db()` for creating BLAST database, and `blast_formatter()` for formatting BLAST output,
the package `blastr` also provide functions for downloading reference genomes from NCBI
(taken from [ncbi.r](https://github.com/J-Moravec/ncbi.r), and a functions to perform a large
number of blast searches at the same time, to improve workflows such as identification of
ortholog genes.

## Installation

Use your favorite of `pak`, `remotes`, `mpd`, or `base::install.packages(repo = NULL)`

Such as using the [mpd](https://github.com/J-Moravec/mpd/) package:
```
mpd::install_github("J-Moravec/blastr")
```

## Dependencies

`blastr` is a binding to [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), so you need to have it installed and available in path.

Follow the installation instructions on the [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) website. You can verify that you have BLAST installed by booting R and typing:

```
Sys.which("blastn")
Sys.which("makeblastdb")
```

If `blastn` is present, other variants (`blastp`, `tblastn`, ...) will most likely be present as well.

## Usage
