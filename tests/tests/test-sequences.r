datadir = "../data"
fna = "GCF_000091225.2.fna.gz"
faa = "GCF_000091225.2.faa.gz"

nuc = read_fasta(file.path(datadir, fna))
prot = read_fasta(file.path(datadir, faa))

TEST_SET("read_fasta correctly reads nucleotide and protein fasta", {
    TEST(inherits(nuc, "sequences"))
    TEST(!is.null(names(nuc)) && names(nuc)[1] |> substring(1, 11) == "NC_003242.2")
    TEST(length(nuc) == 11 && nchar(nuc[[1]]) == 209982)

    TEST(inherits(prot, "sequences"))
    TEST(!is.null(names(prot)) && names(prot)[1] |> substring(1, 14) == "NP_001402090.1")
    TEST(length(prot) == 2122 && nchar(prot[[1]]) == 158)
    })


TEST_SET("is_sequences can correctly guess sequences", {
    n = head(nuc)
    p = head(prot)

    # These should all be detected as sequences
    TEST(is_sequences(nuc) && is_sequences(prot))
    TEST(is_sequences(n) && is_sequences(p))
    TEST(is_sequences(substring(n, 1, 10)) && is_sequences(substring(p, 1, 10)))

    tmp = tempfile()
    TEST(!is_sequences(file.path(datadir, fna))) # full path to existing file
    TEST(!is_sequences(fna)) # no / but extension
    TEST(!is_sequences(tmp)) # no extension, but path separator (either /, \, or \\)

    # Ambiguous case, reported as a sequence
    TEST(is_sequences(basename(tmp))) # ambiguous?
    })
