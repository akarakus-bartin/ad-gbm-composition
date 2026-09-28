# tools/import_downloads.R — move files downloaded from the chat into the project
# Usage:  source(here::here("tools/import_downloads.R"))
#         import_downloads()              # show plan only (safe to paste with other lines)
#         import_downloads(move = TRUE)   # actually move the files
import_downloads <- function(move = FALSE, downloads = "~/Downloads", days = 3) {
  proj <- here::here()
  rules <- data.frame(
    pattern = c("^import_downloads\\.R$",
                "^[0-9]{2}[a-z]?_[A-Za-z0-9_]+\\.R$",
                "^Figure[0-9]+_.*\\.(pdf|png|svg)$",
                "^FigureS[0-9]+_.*\\.(pdf|png|svg)$",
                "_scenarioB_.*\\.(docx|md)$",
                "_NoA\\.(docx|md)$",
                "^facts_.*\\.md$"),
    dest = c("tools", "R", "outputs/figures_scenarioB", "outputs/supplementary",
             "manuscript/drafts", "manuscript/drafts", "manuscript"),
    stringsAsFactors = FALSE)
  f <- list.files(path.expand(downloads), full.names = TRUE)
  f <- f[!file.info(f)$isdir & file.info(f)$mtime > Sys.time() - days * 86400]
  clean <- sub(" \\([0-9]+\\)(\\.[^.]+)$", "\\1", basename(f))
  dest <- vapply(clean, function(x) {
    h <- which(vapply(rules$pattern, grepl, logical(1), x = x))
    if (length(h)) rules$dest[h[1]] else NA_character_ }, character(1))
  keep <- !is.na(dest)
  if (!any(keep)) { message("No project files among recent downloads."); return(invisible(NULL)) }
  plan <- data.frame(from = f[keep], name = clean[keep], dest = dest[keep],
                     mtime = file.info(f[keep])$mtime, stringsAsFactors = FALSE)
  plan <- plan[order(plan$name, plan$mtime, decreasing = TRUE), ]
  plan <- plan[!duplicated(plan$name), ]
  plan$to <- file.path(proj, plan$dest, plan$name)
  print(data.frame(file = plan$name, to = plan$dest,
                   action = ifelse(file.exists(plan$to), "overwrite", "new"),
                   downloaded = format(plan$mtime, "%d.%m %H:%M")), row.names = FALSE)
  if (!move) { message("\nPlan only. Run import_downloads(move = TRUE) to move these files."); return(invisible(plan)) }
  for (i in seq_len(nrow(plan))) {
    dir.create(dirname(plan$to[i]), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(plan$from[i], plan$to[i], overwrite = TRUE, copy.date = TRUE)
    if (ok) file.remove(plan$from[i])
    cat(sprintf("  %s  %s -> %s\n", if (ok) "OK" else "FAILED", plan$name[i], plan$dest[i]))
  }
  left <- f[keep][file.exists(f[keep])]
  if (length(left)) { file.remove(left); cat("  Removed", length(left), "older duplicate(s)\n") }
  invisible(plan)
}

