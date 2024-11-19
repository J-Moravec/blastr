workdir = tempdir()
datadir = "../data"
name = "GCF_000091225.2.fna"

if(!nzchar(Sys.which("makeblastdb"))){
    message("SKIP test-make_blast_db.r -- no makeblastdb binary")
    return(invisible())
    }


TEST_SET("make_blast_db creates db from fasta", {
    file.copy(
        file.path(datadir, paste0(name, ".gz")),
        file.path(workdir, paste0(name, ".gz"))
        )

    input = gunzip(
        file.path(workdir, paste0(name, ".gz")),
        file.path(workdir, name)
        )

    TEST_NOT_ERROR(make_blast_db(input, file.path(workdir, "db", name)))
    dir.remove(file.path(workdir, "db"))
    file.remove(input)
    })


TEST_SET("make_blast_db creates db from gzipped fasta", {
    input = file.path(workdir, paste0(name, ".gz"))
    file.copy(
        file.path(datadir, paste0(name, ".gz")),
        input
        )

    TEST_NOT_ERROR(make_blast_db(input, file.path(workdir, "db", name)))
    dir.remove(file.path(workdir, "db"))
    file.remove(input)
    })


TEST_SET("make_blast_db creates db from sequences", {
    input = read_fasta(file.path(datadir, paste0(name, ".gz")))

    TEST_NOT_ERROR(make_blast_db(input, file.path(workdir, "db", name)))
    })
