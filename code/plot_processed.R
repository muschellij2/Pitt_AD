library(neurobase)
library(scales)
library(dplyr)
setwd(here::here())
ids = list.files(path = "data", recursive = FALSE)
df = tibble::tibble(full_id = ids) %>%
  mutate(
    field_strength = ifelse(grepl("1.5T", full_id), 1.5, 3),
    id = sub("_1.5T", "", full_id),
    outdir = paste0("proc/", id))
i = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(i)) {
  i = 1
}
# for (i in seq(nrow(df))) {
  # print(i)
  idf = df[i,]
  iid = idf$id
  outdir = idf$outdir
  files = c(
    img = "T1_N4_noneck_reduced_winsor_regtoT1.nii.gz",
    brain_mask = "Brain_Mask.nii.gz",
    brain = "T1_N4_noneck_reduced_winsor_regtoT1_brain_N4.nii.gz",
    csf = "CSF.nii.gz",
    gm = "GM.nii.gz",
    wm = "WM.nii.gz",
    tissue = "TISSUES.nii.gz",
    structures = "STRUCTURES.nii.gz",
    cort_thickness = "corticalThickness.nii.gz"
  )
  nm = names(files)
  files = file.path(outdir, files)
  names(files) = nm
  fe = file.exists(unlist(files[c("img", "brain_mask", "tissue", "brain")]))
  files = files[file.exists(files)]
  files = as.list(files)
  alpha_red = "#FF000080"
  
  pdfname = file.path(outdir, "prenormalize_results.pdf")
  
  if (all(fe) & !file.exists(pdfname)) {
    img = readnii(files$img)
    bm = readnii(files$brain_mask)
    pdf(pdfname)
    ortho2(img, bm, col.y = alpha_red)
    tissue = readnii(files$tissue)
    ortho2(img, tissue, col.y = alpha(c("green", "blue", "red"), 0.5))
    brain = readnii(files$brain)
    # cort = readnii(files$cort_thickness)
    # ortho2(brain, cort)
    dev.off()
    outpdf = file.path("results", paste0(iid, "_prenormalize_results.pdf"))
    file.copy(pdfname, outpdf, overwrite = TRUE)
  }
  
# }
  