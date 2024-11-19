datadir = "../data"
name = "GCF_000091225.2.fna.gz"


if(!nzchar(Sys.which("blastn"))){
    message("SKIP test-blast.r -- no blastn binary")
    return(invisible())
    }


TEST_SET("blast produce expected result", {
    subject = read_fasta(file.path(datadir, name))
    query = subject[1] |> substring(1, 100) |> structure(class = "sequences")

    TEST(nrow(blastn(query, subject, outfmt = 6)) == 4)
    })


TEST_SET("blast throws error with wrong type", {
    faa = read_fasta(file.path(datadir, "GCF_000091225.2.faa.gz"))
    TEST_ERROR(blastn(head(faa), head(faa), outfmt = 6, quiet = TRUE))
    })
