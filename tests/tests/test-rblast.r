if(!nzchar(Sys.which("blastp"))){
    message("SKIP test-blast.r -- no blastp binary")
    return(invisible())
    }


# shared data and variables
datadir = "../data"
name = "GCF_000091225.2.faa.gz"

# Notes -- Here are some reasons why your rblast might fail:
# -- non-unique sequences -- then it is essentially a random roll which match is the "best"
# -- sequence is a subset of another -- BLAST algorithm gives a higher score to longer
#                                       subject if it contains a sequence, even if there is
#                                       a subject identical to query
# This might be resolved by a smarter rblast (best n match blast)
sequences = read_fasta(file.path(datadir, name))
sequences = sequences[!duplicated(sequences)] # unique removes names
sequences = head(sequences, 10)


TEST_SET("rblast on itself returns itself as a match", {
    res = rblast(
        head(sequences, 5) |> names(),
        query = sequences,
        subject = sequences,
        type = "blastp"
        )

    TEST(all(res$query == res$subject))
    TEST(all(res$perc_identity == 100))
    TEST(all(res$mismatches == 0))
    })


TEST_SET("rblast reports an error if input is not in query exactly once", {
    x = names(sequences) |> head(5)

    TEST_ERROR(
        rblast(c(x, "foo"), sequences, sequences, type = "blastp"),
        pattern = "Some sequence names were not found: foo"
        )

    TEST_ERROR(
        rblast(x, c(sequences, sequences), sequences, type = "blastp"),
        pattern = "Non unique match for sequence names:"
        )
    })


TEST_SET("rsblast on itself returns itself as a match depending on p", {
    res = rsblast(
        head(sequences, 5) |> names(),
        query = sequences,
        subject = sequences,
        type = "blastp",
        n = 10,
        p = 0.2
        )

    TEST(all(res$query == res$subject))
    TEST(all(res$perc_identity == 100))
    TEST(all(res$mismatches == 0))
    })
