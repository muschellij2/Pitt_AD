library(RNifti)
library(neurobase)
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd(here::here())
ids = list.files(path = "data", recursive = FALSE)
df = tibble::tibble(full_id = ids) %>%
  mutate(
    field_strength = ifelse(grepl("1.5T", full_id), 1.5, 3),
    id = sub("_1.5T", "", full_id),
    proc_dir = paste0("proc/", full_id),
    data_dir = paste0("data/", full_id))
df$original_image = sapply(df$data_dir, list.files, pattern = "mprg|spgr", 
                           full.names = TRUE)

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

for (ifile in names(files)) {
  f = files[ifile]
  df[, ifile] = file.path(df$proc_dir, f)
}

wide = df %>% 
  select(id, field_strength, tissue) %>% 
  spread(field_strength, tissue) %>% 
  rename(weak = "1.5",
         strong = 3)

iid = 1 
results = vector(mode = "list", length = nrow(df))
names(results) = df$full_id
for (iid in seq(nrow(df))) {
  print(iid)
  x = unlist(df[iid, c("brain_mask", "tissue")])
  if (all(file.exists(x))) {
    bm = df$brain_mask[iid]
    bm = readNifti(bm)
    tissue = df$tissue[iid]
    tissue = readNifti(tissue)
    vals = tissue[bm == 1]
    tab = voxres(tissue, "cm")*table(vals)
    tab = as.numeric(tab)
    results[[iid]] = tab
  }
}
xresults = results

results = xresults
results = results[!sapply(results, is.null)]
results = lapply(results, function(x) as.data.frame(t(x)))
mat = dplyr::bind_rows(results, .id = "full_id")
colnames(mat) = c("full_id", "background", "CSF", "GM", "WM")
mat = as.data.frame(mat)
mat$total = rowSums(mat[, c("background", "CSF", "GM", "WM")])

long = right_join(df, mat) %>% 
  select(id, field_strength, background, CSF, GM, WM, total) %>% 
  gather(var, value, background, CSF, GM, WM, total) 
ddf = long %>% 
  mutate(var = paste0(var, "_", field_strength)) %>% 
  select(-field_strength) %>% 
  spread(var, value)

long = long %>% 
  spread(field_strength, value)

ctypes = c( "CSF", "GM", "WM")

pngname = file.path("results", "comparison_of_volumes.png")
png(pngname, height = 3.5, width = 10, units = "in", res = 300)

mod = long %>% 
  filter(var %in% ctypes) %>% 
  ggplot(aes(x = `1.5`, y = `3`)) + geom_point() +
  facet_wrap(~ var, scales = "free") + 
  geom_abline(intercept = 0, slope = 1, col = "black") + 
  geom_smooth(se = FALSE, method = "lm", col = "red") + 
  geom_smooth(se = FALSE) 
print(mod)
dev.off()


long %>% 
  mutate(diff = `1.5` - `3`) %>% 
  group_by(var) %>% 
  summarise(mn = mean(diff),
            sd = sd(diff),
            mn_abs = mean(abs(diff)),
            rmse = sqrt(mean(diff^2)))

res = rep(NA, length(ctypes))
names(res) = ctypes
for (itype in ctypes) {
  cols = paste0(itype, c("_1.5", "_3"))
  res[itype] = cor(ddf[, cols], use = "pairwise.complete.obs")[1,2]
}

res

