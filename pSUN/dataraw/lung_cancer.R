### dataset downloaded at the following link:
# http://archive.ics.uci.edu/dataset/62/lung+cancer
### downloaded 28/08/2023

dataset_lung = read.table("lung-cancer.data", header = FALSE, sep = ",")

y_LungCancer = as.numeric(dataset_lung$V1)
y_LungCancer[y_LungCancer != 1] = 0
X_LungCancer = matrix(nrow = 32, ncol = 54)
X_LungCancer = dataset_lung[, -c(1, 5, 39)]

colnames(X_LungCancer) = c(
  paste("Attribute", c(2:4), sep = ""),
  paste("Attribute", c(6:38), sep = ""),
  paste("Attribute", c(40:57), sep = "")
)

lung_cancer = cbind(y_LungCancer, X_LungCancer)

usethis::use_data(lung_cancer)
