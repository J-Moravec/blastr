# blastr: mapping to NCBI BLAST

`blastr` is a lightweight no-dependency package providing a binding to the [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), a local version of [BLAST](https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi).

Outside of `blast()` and its wrappers `blastn()`, `blastp()`, etc. for blasting,
`make_blast_db()` for creating BLAST database, and `blast_formatter()` for formatting BLAST output (not yet implemented),
the package `blastr` also provide functions for downloading reference genomes from NCBI
(taken from [ncbi.r](https://github.com/J-Moravec/ncbi.r), and a `rblast()` function to perform reciprocal best hit search (RBH), to improve workflows such as identification of
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

Other than R and NCBI BLAST+, no other dependencies are required.

## Usage

Let's first download the data for *Encephalitozoon cuniculi*, a fungal parasite with a very small genome.
Due to their small size, their sequences were used to unit-test `blastr`.

Here we download the protein fasta sequences:

```{r}
ecuni = "GCF_000091225.2"
blastr::ncbi(ecuni, "PROT_FASTA", dir = ".")
```

Now, we want to find the *pho4* gene (like the Vietnamese soup), an important S. cerevisie gene responsible for the phosphate metabolism.

So we copy the part of *pho4* from genbank:

```{r}
# from: https://www.ncbi.nlm.nih.gov/protein/285811932
pho4 = "MGRTTSEGIHGFVDDLEPKSSILDKVGDFITVNTKRHDGREDFNEQNDELNSQENHNSSENGNENENEQD"
names(pho4) = "pho4" # <- needs to be named
```

Now we can run `blastp` to find the matching *pho4* protein (or its part) in E. cuni:

```{r}
ecuni_faa = paste0(ecuni, ".faa.gz")
blastr::blastp(pho4, ecuni_faa, outfmt = 6)
```

The first four columns should look like this:

|query |subject        | perc_identity| alignment_length|
|:-----|:--------------|-------------:|----------------:|
|pho4  |NP_597388.1    |        38.710|               31|
|pho4  |NP_001402271.1 |        57.143|               14|
|pho4  |NP_586162.1    |        27.778|               54|
|pho4  |NP_584778.1    |        54.545|               11|
|pho4  |NP_586498.2    |        34.483|               29|
|pho4  |NP_586199.1    |        27.778|               36|

To me, it doesn't look like we got a great match.
