check_prog = function(prog){
    if(!nzchar(Sys.which(prog)))
        stop("Program \"", prog, "\" is not available.")
    }

