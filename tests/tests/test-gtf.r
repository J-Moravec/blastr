gtf = file.path("../data", "GCF_000091225.2.gtf.gz")

TEST_SET("read_gtf can read files in gtf format", {
    g = read_gtf(gtf)
    TEST(ncol(g) == 9)
    TEST(nrow(g) == 12950) # 1255 lines, 5 of them comments

    # subset for less work
    g = g[g$feature == "gene",]
    TEST(nrow(g) == 2155)

    attr = parse_attributes(g$attribute)
    TEST(ncol(attr) == 9)
    TEST(nrow(attr) == 2155)

    gt = read_gtf(gtf, feature = "gene", attributes = TRUE)
    TEST(ncol(gt) == 17)
    TEST(nrow(gt) == 2155)
    TEST(identical(gt[1:8], g[1:8]))
    # gt, but not attr has rownames from subsetting
    rownames(gt) = NULL
    rownames(attr) = NULL
    TEST(identical(gt[-c(1:8)], attr))
    })
