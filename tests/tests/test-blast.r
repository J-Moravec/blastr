if(!nzchar(Sys.which("blastn"))){
    message("SKIP test-blast.r -- no blastn binary")
    return(invisible())
    }


datadir = "../data"
name = "GCF_000091225.2"
fna = file.path(datadir, paste0(name, ".fna.gz"))
faa = file.path(datadir, paste0(name, ".faa.gz"))
seq = read_fasta(fna)


TEST_SET("blast produce expected result", {
    TEST(nrow(blastn(substring(seq[1], 1, 100), seq, outfmt = 6)) == 4)
    })


TEST_SET("blast throws error with wrong type", {
    faa_seq = read_fasta(faa) |> head()
    TEST_ERROR(blastn(faa_seq, faa_seq, outfmt = 6, quiet = TRUE))
    })


TEST_SET("blast can be fed sequences, files, and gzipped files", {
    seq_head = head(seq, 1)
    fna_head = tempfile(fileext = ".fna")
    write_fasta(seq_head, fna_head)
    fna_head_gz = gzip(fna_head, keep = TRUE)

    # test for query
    TEST_NOT_ERROR(blastn(fna_head, seq_head, quiet = TRUE))
    TEST_NOT_ERROR(blastn(fna_head_gz, seq_head, quiet = TRUE))

    # test for subject
    TEST_NOT_ERROR(blastn(seq_head, fna_head, quiet = TRUE))
    TEST_NOT_ERROR(blastn(seq_head, fna_head_gz, quiet = TRUE))
    file.remove(c(fna_head, fna_head_gz))
    })
