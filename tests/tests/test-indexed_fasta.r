fasta = "../data/seq.fasta"

TEST_SET("indexed_fasta object can be initialized and subsetted", {
    seq = read_indexed_fasta(fasta)

    # Read sequences
    TEST(identical(seq[[1]], structure(c("one" = "AAA"), class = c("sequences", "character"))))
    TEST(identical(unclass(seq[[3]]), c("three" = "TTTTTT")))
    TEST(identical(unclass(seq[[c(2,1)]]), c("two" = "CCC", "one" = "AAA")))

    # Subset
    TEST(identical(seq[1:3], seq))
    TEST(identical(class(seq[3]), "indexed_fasta"))
    TEST(identical(unclass(seq[3][[1]]), c("three" = "TTTTTT")))

    # Subsets with names
    TEST(identical(seq[1], seq["one"]))
    TEST(identical(seq[c(3,1)], seq[c("three", "one")]))
    TEST(identical(seq[[2]], seq[["two"]]))

    # names getter and setter
    TEST(identical(names(seq), c("one", "two", "three")))
    s = setNames(seq, c("foo", "bar", "baz")) # create a copy
    TEST(identical(names(s), c("foo", "bar", "baz")))
    TEST(seq[[1]] == s[[1]])

    # $ operator
    TEST(identical(seq$one, seq[[1]]))

    # Undefined operators should error
    # [<-, [[<-, and $<- don't make sense in this context
    TEST_ERROR({seq$one = "foo"})
    TEST_ERROR({seq[1] = "foo"})
    TEST_ERROR({seq[[1]] = "foo"})
    })
