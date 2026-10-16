# Visual review helper — reconstructed 2026-10-16 (K4 work order
# runs/2026-10-16-test-infra-landing). Referenced by the 2026-09-15 inventory
# but missing from disk; 46 test errors on first baseline run.
#
# Contract (inferred from ~46 call sites): called AFTER assertions with a
# ggplot object + human-readable label; side-effect only; must never make a
# test fail spuriously; must never write into the package directory.
#
# Output dir resolution:
#   1. getOption("seuratextend.review_dir") — set by the main session when
#      collecting plots for human visual review (e.g. runs/<wo>/review-plots/)
#   2. file.path(tempdir(), "seuratextend-review") — default; R CMD check-safe
# Save failure -> warning (not error): this is a review aid, not an assertion.
save_for_review <- function(p, label) {
  dir <- getOption("seuratextend.review_dir",
                   file.path(tempdir(), "seuratextend-review"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(dir, paste0(make.names(label), ".pdf"))
  tryCatch(
    ggplot2::ggsave(file, plot = p, width = 7, height = 5),
    error = function(e) warning("save_for_review failed for '", label, "': ",
                                conditionMessage(e), call. = FALSE)
  )
  invisible(file)
}
