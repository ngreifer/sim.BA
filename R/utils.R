match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg().")
  arg.name <- deparse1(substitute(arg))

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg)) return(choices[1L])

  if (!is.character(arg)) {
    if (several.ok)
      chk::err(sprintf("the argument to `%s` must be `NULL` or a character vector", arg.name))
    else
      chk::err(sprintf("the argument to `%s` must be `NULL` or a string", arg.name))
  }


  if (!several.ok) {
    if (identical(arg, choices)) return(arg[1L])
    if (length(arg) > 1L) {
      chk::err(sprintf("the argument to `%s` must be of length 1", arg.name))
    }
  }
  else if (length(arg) == 0) {
    chk::err(sprintf("the argument to `%s` must be of length >= 1", arg.name))
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    chk::err(sprintf("the argument to `%s` should be %s%s",
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2)))

  i <- i[i > 0L]

  choices[i]
}

word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  word.list <- add_quotes(word.list, quotes)

  if (L == 0) {
    out <- ""
    attr(out, "plural") <- FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") <- FALSE
    }
    else {
      and.or <- match_arg(and.or, c("and", "or"))
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or, " "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L - 1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or, " "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") <- TRUE
    }

  }

  out
}

add_quotes <- function(x, quotes = 2L) {
  if (!isFALSE(quotes)) {
    if (isTRUE(quotes)) quotes <- 2

    if (chk::vld_string(quotes)) x <- paste0(quotes, x, quotes)
    else if (chk::vld_whole_number(quotes)) {
      if (as.integer(quotes) == 0) return(x)
      else if (as.integer(quotes) == 1) x <- paste0("\'", x, "\'")
      else if (as.integer(quotes) == 2) x <- paste0("\"", x, "\"")
      else stop("`quotes` must be boolean, 1, 2, or a string.")
    }
    else {
      stop("'quotes' must be boolean, 1, 2, or a string.")
    }
  }
  x
}

# Capitalize first letter of string
firstup <- function(x) {
  `substr<-`(x, 1, 1, toupper(substr(x, 1, 1)))
}