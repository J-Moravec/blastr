mock = function(x, val, env){
    orig = get(x, envir = env)
    unlockBinding(x, env)
    assign(x, val, env)
    lockBinding(x, env)

    invisible(orig)
    }


with_mock = function(x, val, env, ...){
    orig = mock(x, val, env)
    on.exit(
        mock(x, orig, env), add = TRUE
        )

    eval(...)
    }


with_mock_pkg = function(x, val, pkg, ...){
    ns = getNamespace(pkg)
    orig = mock(x, val, ns)
    on.exit(mock(x, orig, ns), add = TRUE)

    env = pkg_env(pkg)
    if(!is.null(env)){
        mock(x, val, env)
        on.exit(mock(x, orig, env), add = TRUE)
        }

    eval(...)
    }


pkg_env = function(pkg){
    if(!is.character(pkg))
        stop("Provide package name")

    nm = paste0("package:", pkg[1])
    m = match(nm, search(), nomatch = 0)

    if(m) pos.to.env(m) else NULL
    }
