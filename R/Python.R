load_condaenv <- function(conda_env) {
  library(reticulate)

  # # Check if Conda is installed
  #
  # if (!conda_exists()) {
  #   # Install Miniconda if Conda is not installed
  #   message("Conda is not installed. Installing Miniconda...")
  #   install_miniconda()
  # }

  # Check if the specified Conda environment exists
  envs <- conda_list()
  if (!conda_env %in% envs$name) {
    if (conda_env == "seuratextend") {
      # Ask the user for confirmation to create the 'seuratextend' environment
      message(paste(
        "Conda environment 'seuratextend' does not exist.\n",
        "Would you like to create it with packages Scanpy, CellRank, scVelo, Palantir?",
        "This environment will require approximately 4GB of space. (yes/no): "
      ))
      response <- tolower(readline(prompt = ""))

      if (response %in% c("yes","y")) {
        create_condaenv_seuratextend()
        use_condaenv(conda_env, required = TRUE)
        return(TRUE)
      } else {
        message("Aborted by the user.")
        return(FALSE)
      }
    } else {
      stop(paste(
        glue::glue("The Conda environment '{conda_env}' does not exist."),
        "Please use 'seuratextend' or create the environment manually."
      ))
    }
  }

  # Load the specified Conda environment
  use_condaenv(conda_env, required = TRUE)
  py_run_string("")
  return(TRUE)
}

# conda_exists <- function() {
#   # Try to find the path to the conda executable
#   conda_path <- Sys.which("conda")
#
#   # Check if the conda executable exists
#   !is.na(conda_path) && file.exists(conda_path)
# }

#' @title Create the 'seuratextend' Conda Environment
#' @description Initializes a Conda environment named 'seuratextend' which includes all necessary Python packages for running `scVelo`, `Palantir`, and `CellRank` via `SeuratExtend`. This function is ideal for users setting up their analytical environment for the first time or those who need to ensure they have a compatible environment setup.
#' @return Does not return any value; it creates and configures a new Conda environment.
#' @details When you run functions from `SeuratExtend` that require Python integrations, such as `scVelo`, `Palantir`, or `CellRank`, you might need a specific set of Python packages installed. This function facilitates the setup of a Conda environment specifically tailored for these tools, managed through `reticulate`. The environment setup is automatic and checks your operating system to configure appropriately. It is supported on Windows, macOS, and Linux (Ubuntu 20.04), ensuring broad compatibility and ease of setup.
#' @examples
#' \dontrun{
#' # To manually create the 'seuratextend' Conda environment, execute:
#' create_condaenv_seuratextend()
#' }
#' @rdname create_condaenv_seuratextend
#' @export

create_condaenv_seuratextend <- function() {
  library(reticulate)

  # Define the package name and the environment name
  env_name <- "seuratextend"

  # Check the operating system
  os_type <- Sys.info()["sysname"]

  # Define the YAML file based on the operating system
  yaml_file <- system.file(
    "extdata",
    switch (os_type,
            "environment-linux.yml",
            Windows = "environment-windows.yml",
            Darwin = "environment-mac.yml"
            ),
    package = "SeuratExtend",
    mustWork = TRUE)

  # Create the Conda environment using the specified YAML file
  message(glue::glue("Creating Conda environment '{env_name}' using '{yaml_file}'..."))
  conda_create(envname = env_name, environment = yaml_file)

  message(glue::glue("Conda environment '{env_name}' has been created successfully."))
}

# check if OS is Windows
is_windows <- function() {
  os_type <- Sys.info()["sysname"]
  return(os_type == "Windows")
}

normalizePathWindows <- function(path) {
  if(is_windows()) path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  return(path)
}

create_temp_dir <- function() {
  # Create a temporary directory with a timestamp
  tmp_dir <- tempdir()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  subfolder <- file.path(tmp_dir, paste0("PalantirPs_", timestamp))
  subfolder <- normalizePathWindows(subfolder)
  dir.create(subfolder)
  return(subfolder)
}

# Helper function to convert R vector to Python list or tuple string
r_vector_to_py <- function(vec, type = "list") {
  if (is.null(vec)) return(NULL)
  start_char <- ifelse(type == "list", "[", "(")
  end_char <- ifelse(type == "list", "]", ")")
  # Check if the vector elements are numeric
  if (is.numeric(vec)) {
    paste0(start_char, paste(vec, collapse = ", "), end_char)
  } else {
    paste0(start_char, paste(shQuote(vec, type = "sh"), collapse = ", "), end_char)
  }
}
