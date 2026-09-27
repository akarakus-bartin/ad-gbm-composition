# =============================================================================
# 07c_allen_reference_subset.R — select reference cells for MuSiC (D13-D15, D20-D21)
# Uses metadata only (no expression, no bulk outcome).
# =============================================================================
d  <- here::here("data/raw/allen_m1")
md <- data.table::fread(file.path(d, "metadata.csv"))
map <- c(Astro = "ast", `Micro-PVM` = "mic", Oligo = "oli", OPC = "opc", Endo = "end")
md$ref_class <- ifelse(md$class_label %in% c("Glutamatergic", "GABAergic"), "neu",
                       unname(map[md$subclass_label]))
md <- md[!is.na(md$ref_class), ]
set.seed(PARAMS$seed)
sel <- do.call(rbind, lapply(split(md, list(md$ref_class, md$external_donor_name_label)), function(x)
  if (nrow(x) > 500) x[sample(nrow(x), 500), ] else x))
sel <- sel[, c("sample_name", "ref_class", "external_donor_name_label")]
print(table(sel$ref_class, sel$external_donor_name_label))
data.table::fwrite(sel, file.path(d, "reference_cells_selected.csv"))
writeLines(sel$sample_name, file.path(d, "keep_ids.txt"))
message("Selected ", nrow(sel), " cells; wrote keep_ids.txt")

