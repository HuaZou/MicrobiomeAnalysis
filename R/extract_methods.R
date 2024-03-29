#' @title Extract `marker_table` object
#'
#' @description Operators acting on `marker_table` to extract parts.
#'
#' @name [
#' @aliases [,marker_table,ANY,ANY,ANY-method
#' @inheritParams base::Extract
#' @param ... see [`base::Extract()`].
#' @param drop ignored now.
#' @return a `marker_table` object.
#' @export
#' @rdname extract-methods
#' @seealso [`base::Extract()`]
#'
setMethod("[", "marker_table", function(x, i, j, ...) {
    newx <- marker_table(data.frame(x)[i, j, drop = FALSE])
    newx
})
