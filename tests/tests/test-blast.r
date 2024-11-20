datadir = "../data"
name = "GCF_000091225.2"
fna = file.path(datadir, paste0(name, ".fna.gz"))
faa = file.path(datadir, paste0(name, ".faa.gz"))

if(!nzchar(Sys.which("blastn"))){
    message("SKIP test-blast.r -- no blastn binary")
    return(invisible())
    }


TEST_SET("blast produce expected result", {
    seq = read_fasta(fna) 
    TEST(nrow(blastn(substring(seq[1], 1, 100), seq, outfmt = 6)) == 4)
    })


TEST_SET("blast throws error with wrong type", {
    seq = read_fasta(faa)
    TEST_ERROR(blastn(head(seq), head(seq), outfmt = 6, quiet = TRUE))
    })


TEST_SET("blast can be fed sequences, files, and gzipped files", {
    seq = read_fasta(fna)
    seq_head = head(seq)
    # gunzipped fna, gunziping is done on a copy just in case
    gufna = temp.copy(fna, fileext = ".fna.gz") |> gunzip()

    # test for query
    TEST_NOT_ERROR(blastn(seq_head, seq_head, quiet = TRUE))
    TEST_NOT_ERROR(blastn(fna, seq_head, quiet = TRUE))
    TEST_NOT_ERROR(blastn(gufna, seq_head, quiet = TRUE))

    # test for subject
    TEST_NOT_ERROR(blastn(seq_head, fna, quiet = TRUE))
    TEST_NOT_ERROR(blastn(seq_head, gufna, quiet = TRUE))
    file.remove(gufna)
    })
