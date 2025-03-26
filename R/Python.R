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
#' @param force Logical; if TRUE, removes any existing 'seuratextend' environment before creating a new one. Default is FALSE. Use this parameter to force recreation of the environment if a previous installation was interrupted or incomplete.
#' @return Does not return any value; it creates and configures a new Conda environment.
#' @details When you run functions from `SeuratExtend` that require Python integrations, such as `scVelo`, `Palantir`, or `CellRank`, you might need a specific set of Python packages installed. This function facilitates the setup of a Conda environment specifically tailored for these tools, managed through `reticulate`.
#'
#' The environment setup is automatic and checks your operating system to configure appropriately. It is currently supported and tested on:
#' * Windows (Intel/AMD processors)
#' * macOS (both Intel-based and Apple Silicon/M-series processors)
#' * Linux (Ubuntu 20.04)
#'
#' The function automatically detects Apple Silicon (M1/M2/M3/M4) Macs and uses the appropriate configuration.
#'
#' If force=TRUE is used, any existing 'seuratextend' environment will be removed before creating a new one. This can be useful when:
#' * A previous installation was interrupted
#' * Required dependencies were not properly installed
#' * You need to reset the environment to a clean state
#' * You encounter persistent environment-related errors
#' @examples
#' \dontrun{
#' # Standard installation
#' create_condaenv_seuratextend()
#'
#' # Force recreation of the environment
#' create_condaenv_seuratextend(force = TRUE)
#' }
#' @rdname create_condaenv_seuratextend
#' @export

create_condaenv_seuratextend <- function(force = FALSE) {
  library(reticulate)

  # Check if git is installed
  git_path <- Sys.which("git")
  if (git_path == "") {
    warning(paste(
      "Git appears to not be installed on your system.",
      "You may need to install Git to properly set up the conda environment.",
      "Please install Git based on your operating system and restart R if you encounter installation issues.",
      "\nInstallation instructions:",
      "\n- Windows: Download from https://gitforwindows.org",
      "\n- Mac: Run 'brew install git' or 'xcode-select --install'",
      "\n- Linux (Ubuntu/Debian): Run 'sudo apt-get install git'",
      "\n- Linux (Fedora): Run 'sudo dnf install git'",
      "\n- Linux (CentOS/RHEL): Run 'sudo yum install git'"
    ))
  }

  # Define the package name and environment name
  env_name <- "seuratextend"

  # Check if environment already exists
  envs <- conda_list()
  if (env_name %in% envs$name) {
    if(!force) {
      message(paste(
        "Conda environment 'seuratextend' already exists.\n",
        "Options:\n",
        "1) Skip creation (enter: skip)\n",
        "2) Force recreate environment (enter: force)\n",
        "Please choose (skip/force): "
      ))
      response <- tolower(readline(prompt = ""))

      if (response == "skip") {
        message("Skipping environment creation.")
        return(invisible())
      } else if (response != "force") {
        message("Invalid response. Skipping environment creation.")
        return(invisible())
      }
    }

    # If force selected, remove existing environment
    message("Removing existing environment...")
    conda_remove(env_name)
  }

  # Check operating system and get appropriate YAML file
  os_type <- Sys.info()["sysname"]

  # Check if on Apple Silicon Mac
  is_apple_silicon <- FALSE
  if (os_type == "Darwin") {
    # Check architecture - Apple Silicon Macs use arm64 or aarch64
    system_arch <- system("uname -m", intern = TRUE)
    is_apple_silicon <- tolower(system_arch) %in% c("arm64", "aarch64")
  }

  # Select appropriate YAML file based on OS and architecture
  if (os_type == "Darwin" && is_apple_silicon) {
    yaml_file_name <- "environment-mac-silicon.yml"
  } else {
    yaml_file_name <- switch(os_type,
                             "environment-linux.yml",  # Default for Linux
                             Windows = "environment-windows.yml",
                             Darwin = "environment-mac.yml"
    )
  }

  yaml_file <- system.file(
    "extdata",
    yaml_file_name,
    package = "SeuratExtend",
    mustWork = TRUE
  )

  # Create the conda environment
  message(glue::glue("Creating conda environment '{env_name}' using '{yaml_file}'..."))
  tryCatch({
    conda_create(envname = env_name, environment = yaml_file)
    message(glue::glue("Conda environment '{env_name}' has been created successfully."))
  }, error = function(e) {
    message(paste(
      "Error creating conda environment.",
      "If you previously had a partially installed environment,",
      "try running this function again with force=TRUE",
      "\nError details:", e$message
    ))
  })
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

#' Activate Python Environment for SeuratExtend
#'
#' This function activates a Python environment (default: "seuratextend")
#'
#' @param conda_env Character string specifying the conda environment name. Default is "seuratextend".
#' @param verbose Logical indicating whether to print status messages. Default is TRUE.
#' @param packages Character vector of Python packages to import.
#'
#' @return Invisible TRUE if successful
#' @export
#'
#' @examples
#' \dontrun{
#' activate_python()
#' activate_python(conda_env = "my_custom_env", verbose = TRUE)
#' }
activate_python <- function(conda_env = "seuratextend",
                            verbose = TRUE,
                            packages = c("os")) {

  load_condaenv(conda_env = conda_env)

  tryCatch({
    # Import basic packages
    for (pkg in packages) {
      if (verbose) message(paste0("Importing ", pkg, "..."))
      reticulate::py_run_string(paste0("import ", pkg))
    }
    invisible(TRUE)
  },
  error = function(e) {
    stop(paste0("Error activating Python environment: ", e$message))
  })
}
