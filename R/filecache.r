#' Call function with on-disk caching
#'
#' This will perform on-disk of function `f`.
#' The function will be executed only if its output is not present in the cache.
#' Otherwise, the content of the cache is returned.
#'
#' This will perform on-disk cachingg of function `f`.
#' First, the parameters are expanded and hashed using the `digest::digest` package.
#' Then the cache,stored as a RDS file is read and examined.
#' If it contains the hash, the content of the hash is returned.
#' Otherwise, the function is evaluated, the content is saved cached, and returned.
#'
#' Note that the function itself is _not_ hashed. This is to avoid potential problems
#' with different versions of a function that provides identical output.
#' Use different file caches for different functions.
#'
#' @param cache A path to cache file, a `.rds` file containing an environment, or an empty path.
#' @param f A function for which the result will be cached, functions themselves are not hashed,
#' see details and examples.
#' @param ... Arguments to be passed to `f`. Arguments are forced to evaluate any promises and
#' then the argument list is hashed.
#' @return A value resulting from calling `f` with arguments `...`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' myfun = function(x) Sys.sleep(x)
#' filecache("cache.rds", myfun, 3)
#' filecache("cache.rds", myfun, 3) # shouldn't take any time
#' filecache("cache.rds", myfun, 2) # stored in the same cache
#' filecache("cache.rds", myfun, 2)
#'
#' # dont do this
#' filecache("cache.rds", sum, 3) # returns incorrect result
#'                                # since hashed 3 is already in cache
#' }
#'
#' @seealso
#' The [digest::digest()] function used to hash arguments.
#' [base::saveRDS()] and [base::readRDS()] for saving and reading the cache.
#' Packages such as `memoize` or `R.cache` implement much more comprehensive caching.
filecache = function(cache, f, ...){
    if(!requireNamespace("digest", quietly = TRUE))
        stop("Please install the \"digest\" package.")

    # this should force/expand ..., don't substitute!
    args = list(...)
    hash = digest::digest(args)

    memory = new.env(parent = emptyenv())
    if(!dir.exists(dirname(cache)))
        dir.create(dirname(cache))

    if(file.exists(cache))
        memory = readRDS(cache)

    if(utils::hasName(memory, hash))
        return(memory[[hash]])

    res = do.call(f, args)
    memory[[hash]] = res
    saveRDS(memory, cache)

    res
    }

