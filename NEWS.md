# blastr 0.0.3

* New `is_gzip()` function to test if a file is gzipped
* New subset `[` function for the S3 class `sequences` that preserves class
* New `print` function for the S3 class `sequences`
* New class `indexed_fasta` that tries to hard really hard that it is just
  a table and a file path with `[`, `[[`, `$` methods.
* New function `read_indexed_fasta` to read fasta using a `fai` index.
  The actual reading is done using the `[[` method, i.e.,
  `seq[[c(5:8)]]` will read the 5, 6, 7, and 8th sequences from the
  fasta file.

# blastr 0.0.2

* Added `qgene` and `sgene` to the `orthologs()` function to specify
  the column where gene identifier is stored. While the default
  `gene_id` should be accurate for most annotations, sometimes
  `gene` needs to be used instead.

# blastr 0.0.1

* Initial version.
