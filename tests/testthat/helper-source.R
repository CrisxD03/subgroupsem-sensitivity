for (f in list.files(testthat::test_path("..", "..", "R"), pattern = "\\.R$", full.names = TRUE))
  source(f)

python_available <- function() {
  ok <- tryCatch(init_python(), error = function(e) FALSE)
  isTRUE(ok)
}
