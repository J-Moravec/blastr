TEST_SET("strfsplit works with corner cases", {

    TEST(identical(strfsplit("foo bar", " ", fixed = TRUE), list("foo", "bar")))
    TEST(identical(strfsplit("foo", " ", fixed = TRUE), list("foo", NA_character_)))
    TEST(identical(strfsplit(NA, " ", fixed = TRUE), list(NA_character_, NA_character_)))

    # test longer character
    TEST(identical(strfsplit("foo bar baz", " bar ", fixed = TRUE), list("foo", "baz")))

    # test vector
    TEST(identical(
        strfsplit(c("foo bar baz", NA, "foo bar"), " ?bar ?"),
        list(c("foo", NA, "foo"), c("baz", NA, ""))
        ))
    })

