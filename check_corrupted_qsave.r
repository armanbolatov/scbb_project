library(qs)

check_qsave_files <- function(folder_path) {
  files <- list.files(folder_path, pattern = "\\.qsave$", full.names = TRUE)
  corrupted_files <- c()

  for (file in files) {
    tryCatch({
      qread(file) 
    }, error = function(e) {
      message(paste("Corrupted file:", file))
      corrupted_files <<- c(corrupted_files, file)
    })
  }

  if (length(corrupted_files) == 0) {
    message("All .qsave files are valid!")
  } else {
    message("The following .qsave files are corrupted:")
    print(corrupted_files)
  }
}

check_qsave_files("/l/users/darya.taratynova/scbb_project/downloaded")
