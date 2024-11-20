dir.size = function(x, recursive = TRUE){
    stopifnot(is.character(x))
    list.files(x, full.names = TRUE, recursive = recursive) |>
        sapply(file.size) |> sum()
    }


dir.remove = function(x, force = FALSE){
    unlink(x, recursive = TRUE, force = force)
    }


temp.copy = function(x, ...){
    copy = tempfile(...)
    file.copy(x, copy)
    copy
    }
