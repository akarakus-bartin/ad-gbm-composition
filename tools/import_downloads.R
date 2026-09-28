# =============================================================================
# tools/import_downloads.R — move files downloaded from the chat into the project
# Usage (R console, project open):
#   source(here::here("tools/import_downloads.R"))
#   import_downloads()                 # asks for confirmation before moving
#   import_downloads(days = 7)         # look further back
#   import_downloads(dry_run = TRUE)   # only show what would be moved
# Only known file-name patterns are moved; everything else in Downloads is ignored.
# =============================================================================
import_downloads <- function(downloads = "~/Downloads", days = 3, dry_run = FALSE) {
  proj <- here::here()
  rules <- data.frame(
    pattern = c("^import_downloads\\.R$",
                "^[0-9]{2}[a-z]?_[A-Za-z0-9_]+\\.R$",
                "^Figure[0-9]+_.*\\.(pdf|png|svg)$",
                "_scenarioB_.*\\.(docx|md)$",
                "^facts_.*\\.md$"),
    dest    = c("tools", "R", "outputs/figures_scenarioB", "manuscript/drafts", "manuscript"),
    stringsAsFactors = FALSE)

  dl <- path.expand(downloads)
  f  <- list.files(dl, full.names = TRUE)
  f  <- f[!file.info(f)$isdir & file.info(f)$mtime > Sys.time() - days * 86400]
  if (!length(f)) { message("No recent files in ", dl); return(invisible(NULL)) }

  clean <- sub(" \\([0-9]+\\)(\\.[^.]+)$", "\\1", basename(f))
  dest  <- vapply(clean, function(x) {
    hit <- which(vapply(rules$pattern, grepl, logical(1), x = x))
    if (length(hit)) rules$dest[hit[1]] else NA_character_
  }, character(1))
  keep <- !is.na(dest)
  if (!any(keep)) { message("No project files found among recent downloads."); return(invisible(NULL)) }

  plan <- data.frame(from = f[keep], name = clean[keep], dest = dest[keep],
                     mtime = file.info(f[keep])$mtime, stringsAsFactors = FALSE)
  plan <- plan[order(plan$name, plan$mtime, decreasing = TRUE), ]
  plan <- plan[!duplicated(plan$name), ]
  plan$to <- file.path(proj, plan$dest, plan$name)
  plan$action <- ifelse(file.exists(plan$to), "overwrite", "new")

  cat("\nPlanned moves (newest copy of each file):\n")
  print(data.frame(file = plan$name, to = plan$dest, action = plan$action,
                   downloaded = format(plan$mtime, "%d.%m %H:%M")), row.names = FALSE)
  if (dry_run) return(invisible(plan))
  if (interactive() && tolower(readline("\nMove these files? [y/N]: ")) != "y") {
    message("Cancelled."); return(invisible(plan))
  }

  for (i in seq_len(nrow(plan))) {
    dir.create(dirname(plan$to[i]), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(plan$from[i], plan$to[i], overwrite = TRUE, copy.date = TRUE)
    if (ok) file.remove(plan$from[i])
    cat(sprintf("  %s  %s -> %s\n", if (ok) "OK" else "FAILED", plan$name[i], plan$dest[i]))
  }
  dups <- f[keep][!f[keep] %in% plan$from & file.exists(f[keep])]
  if (length(dups)) { file.remove(dups); cat("  Removed", length(dups), "older duplicate(s) from Downloads\n") }
  invisible(plan)
}

