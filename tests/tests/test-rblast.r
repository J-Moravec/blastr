if(!nzchar(Sys.which("blastn"))){
    message("SKIP test-blast.r -- no blastn binary")
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

    TEST(all(res$query == res$ortholog))
    TEST(all(res$perc_identity == 100))
    TEST(all(res$mismatches == 0))
    })


TEST_SET("rblast reports an error if input is not in query", {
    x = names(seq) |> head()
    x = c(x, "foo")

    TEST_ERROR(rblast(x, sequences, sequences, type = "blastp"))
    })
