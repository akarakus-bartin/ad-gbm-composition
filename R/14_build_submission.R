# =============================================================================
# 14_build_submission.R — assemble the journal submission package in one folder
# Output: submission/NoA_<date>/  (numbered files in upload order) + MANIFEST.txt
# Run:    source(here::here("R/14_build_submission.R"))
# =============================================================================
proj <- here::here()
out  <- file.path(proj, "submission", paste0("NoA_", format(Sys.Date(), "%Y-%m-%d")))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fig <- file.path(proj, "outputs/figures_scenarioB"); sup <- file.path(proj, "outputs/supplementary")
drafts <- file.path(proj, "manuscript/drafts")

items <- list(
  c("00_Cover_letter.docx",                       file.path(drafts, "cover_letter_NoA.docx")),
  c("01_Manuscript.docx",                          file.path(drafts, "manuscript_NoA.docx")),
  c("02_Highlights.docx",                          file.path(drafts, "highlights_NoA.docx")),
  c("03_Figure1.pdf",                              file.path(fig, "Figure1_study_design.pdf")),
  c("04_Figure2.pdf",                              file.path(fig, "Figure2_H1.pdf")),
  c("05_Figure3.pdf",                              file.path(fig, "Figure3_H2_RRHO2.pdf")),
  c("06_Figure4.pdf",                              file.path(fig, "Figure4_negative_control.pdf")),
  c("07_Figure5.pdf",                              file.path(fig, "Figure5_DD_sensitivity.pdf")),
  c("10_Supplementary_Tables_S1-S4.xlsx",          file.path(sup, "Supplementary_Tables_S1-S4.xlsx")),
  c("11_Supplementary_Figure_S1.pdf",              file.path(sup, "FigureS1_H1_replication.pdf")),
  c("12_Graphical_abstract.tiff",                 file.path(drafts, "graphical_abstract_NoA.tiff")))

status <- data.frame(file = character(), source = character(), ok = logical(), stringsAsFactors = FALSE)
for (it in items) {
  ok <- file.exists(it[2]) && file.copy(it[2], file.path(out, it[1]), overwrite = TRUE)
  status[nrow(status) + 1, ] <- list(it[1], sub(paste0(proj, "/"), "", it[2], fixed = TRUE), ok)
}

# Supplementary File 1: the analysis plan exactly as committed before analysis (3fcf8d1)
sf1 <- file.path(out, "08_Supplementary_File_1_Analysis_Plan.docx")
rc <- system2("git", c("-C", shQuote(proj), "show", "3fcf8d1:docs/Analysis_Plan_v1.docx"), stdout = sf1)
status[nrow(status) + 1, ] <- list(basename(sf1), "git show 3fcf8d1:docs/Analysis_Plan_v1.docx",
                                   rc == 0 && file.exists(sf1) && file.size(sf1) > 0)

# Supplementary File 2: deviation register, converted from Markdown to Word
sf2 <- file.path(out, "09_Supplementary_File_2_Deviation_Register.docx")
ok2 <- FALSE
if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
  ok2 <- tryCatch({ rmarkdown::pandoc_convert(file.path(proj, "DEVIATIONS.md"), to = "docx", output = sf2); file.exists(sf2) },
                  error = function(e) { message("pandoc error: ", conditionMessage(e)); FALSE })
}
status[nrow(status) + 1, ] <- list(basename(sf2), "DEVIATIONS.md (pandoc -> docx)", ok2)
ann <- file.path(drafts, "supplementary_file2_NoA.docx")
if (file.exists(ann)) { file.copy(ann, sf2, overwrite = TRUE); status[nrow(status), ] <- list(basename(sf2), "manuscript/drafts/supplementary_file2_NoA.docx (annotated guide)", TRUE) }

# [CHECK] fields still open in the manuscript
chk <- NA
m <- file.path(out, "01_Manuscript.docx")
if (file.exists(m)) {
  td <- tempfile(); utils::unzip(m, files = "word/document.xml", exdir = td)
  x <- paste(readLines(file.path(td, "word/document.xml"), warn = FALSE), collapse = "")
  chk <- lengths(regmatches(x, gregexpr("[CHECK", x, fixed = TRUE)))
}

status <- status[order(status$file), ]
man <- c(sprintf("Submission package — built %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
         sprintf("git HEAD: %s", paste(system2("git", c("-C", shQuote(proj), "rev-parse", "--short", "HEAD"), stdout = TRUE), collapse = "")),
         sprintf("Open [CHECK] fields in manuscript: %s", chk), "",
         sprintf("%-50s %-8s %s", "file", "md5", "source"))
for (i in seq_len(nrow(status))) {
  p <- file.path(out, status$file[i])
  md5 <- if (status$ok[i]) substr(unname(tools::md5sum(p)), 1, 8) else "MISSING"
  man <- c(man, sprintf("%-50s %-8s %s", status$file[i], md5, status$source[i]))
}
man <- c(man, "", "Not produced by this script (to be prepared separately):",
         "  - Declaration of competing interest (.docx from Elsevier declarations tool) - required")
writeLines(man, file.path(out, "MANIFEST.txt")); cat(man, sep = "\n")
if (any(!status$ok)) message("\n!! Missing items above — see 'source' column.") else message("\nAll files present in: ", out)
