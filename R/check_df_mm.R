#' Check the degrees of freedom required by the moment-matching method
#'
#' The moment-matching approximations of \code{\link{ptdiff}} and of two
#' continuous endpoints require the fourth moments of both t distributions to
#' exist.
#'
#' @param df_t,df_c Degrees of freedom of the distributions of the treatment
#'   and control means.
#' @param method One of \code{"NI"}, \code{"MM"} or \code{"MC"}.
#' @param alternative The methods suggested in the error message.
#'
#' @return \code{NULL}, invisibly.
#' @noRd
check_df_mm <- function(df_t, df_c, method,
                        alternative = "method = \"NI\" or \"MC\"") {
  if (method == "MM" && (df_t <= 4 || df_c <= 4)) {
    stop(sprintf(paste("method = \"MM\" requires the distributions of both",
                       "group means to have more than 4 degrees of freedom",
                       "(here %s and %s). Use %s."),
                 format(df_t), format(df_c), alternative), call. = FALSE)
  }
  invisible(NULL)
}
