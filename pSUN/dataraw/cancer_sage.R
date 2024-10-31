### dataset downloaded at the following link:
# https://webusers.i3s.unice.fr/~pasquier/web/?Research_Activities___Dataset_Downloads___Cancer_SAGE
### downloaded 01/08/2023

dataset_gene = read.csv("dataset_74-516.csv", header = TRUE, sep = "")
y_CancerSAGE = c(
  0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1,
  0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1
)
X_CancerSAGE = dataset_gene[, -1]
cancer_sage = cbind(y_CancerSAGE, X_CancerSAGE)

dimnames(cancer_sage)[[2]][1] = "positivityDiagnosis"

usethis::use_data(cancer_sage)