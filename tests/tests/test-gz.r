workdir = tempdir()

sample_content = function(file){
    content = "Hello world!"

    writeLines(content, file)

    content
    }


is_gzipped = function(x){
    first_two_bytes = readBin(x, "raw", 2)

    magic_number = charToRaw("\x1F\x8b")

    identical(first_two_bytes, magic_number)
    }


TEST_SET("gzip can compress files", {
    input = file.path(workdir, "input.txt")
    output = file.path(workdir, "output.txt.gz")
    content = sample_content(input)

    gzip(input, output, keep = TRUE)
    TEST(is_gzipped(output))
    TEST(identical(readLines(output), content))
    TEST(file.exists(input))
    file.remove(output)

    content = sample_content(input)

    output = gzip(input)
    TEST(is_gzipped(output))
    TEST(identical(readLines(output), content))
    TEST(!file.exists(input))
    file.remove(output)
    })


TEST_SET("gunzip can decompress files", {
    sample = file.path(workdir, "sample.txt")
    content = sample_content(sample)
    gzip(sample)

    sample_gz = paste0(sample, ".gz")
    TEST(file.exists(sample_gz))

    gunzip(sample_gz)
    TEST(identical(readLines(sample), content))
    TEST(!is_gzipped(sample))
    file.remove(sample)
    })
