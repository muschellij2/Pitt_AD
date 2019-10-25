library(neurobase)
library(extrantsr)
library(smri.process)
setwd(here::here())

ids = list.files(path = "data", recursive = FALSE)
# for (iid in ids) {
i = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(i)) {
  i = 1
}
iid = ids[i]
outdir = paste0("proc/", iid)
x = list.files(path = paste0("data/", iid),
               pattern = "mprg|spgr", full.names = TRUE)
stopifnot(length(x) == 1)
names(x) = c("T1")
x = as.list(x)
tissue = c(
  CSF = "CSF.nii.gz",
  GM = "GM.nii.gz",
  WM = "WM.nii.gz",
  TISSUES = "TISSUES.nii.gz",
  STRUCTURES = "STRUCTURES.nii.gz"
)
n = names(tissue)
tissue = file.path(outdir, tissue)
names(tissue) = n
files = c(tissue, sub(".nii", "_resampled.nii", tissue)) 
tissue = as.list(tissue)
if (!all(file.exists(files))) {
  processed = smri_prenormalize(
    x, 
    outdir = outdir,
    num_templates = 35,
    brain_extraction_method = "robex")
  all_resampled = seg_normalize(
    prenormalize = processed, template = "none", verbose = 2)
  
  # tissue = all_resampled$native$tissue
}
outfile = file.path(outdir, "corticalThickness.nii.gz")
if (!file.exists(outfile)) {
  res = cort_thickness(seg = tissue$TISSUES,
                       gray = tissue$GM,
                       white = tissue$WM,
                       verbose = TRUE)
  
  write_nifti(res, outfile)
}


