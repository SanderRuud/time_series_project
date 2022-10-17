# Play with the data (start working on the project)

library(R.matlab)
setwd("/home/ajo/gitRepos/time_series_project")
data <- readMat("Mouse28-140313_BS0150_HMMready.mat")
dim(data)
length(data)
str(data)

angdata <- data$resampledAwakeHeadAngleData
dim(angdata)
str(angdata)

celldata <- data$celldata
dim(celldata)
str(celldata)
# Not sure what exactly this data contains / what it informs us about?!

# Make train and test data for our model. Done "sequentially" in time to avoid data leakage.
N <- length(celldata[1,])
train_indices <- as.integer(0.8*N) # The first 80% of observations are for training. 
cell_train <- t(celldata[,1:train_indices])
cell_test <- t(celldata[,train_indices:N])
ang_train <- angdata[1:train_indices]
ang_test <- angdata[train_indices:N]

# Fitting a HMM can e.g. be done with https://cran.r-project.org/web/packages/depmixS4/index.html
library(depmixS4)
?depmix

# Hyperparameters and priors.
nstates <- 2 # Fixed number of states. 
pi <- rep(1, nstates)/nstates # Uniform prior state probabilities. 

# Try to fit a HMM. 
cell_train_df <- as.data.frame(cell_train)
colnames(cell_train_df) <- c("x", "y")
model <- depmix(response = y ~ x, family = poisson(link = "log"), data = cell_train_df[1:100,c(1,2)], 
                nstates = nstates) 
# Vi tester med lite data først. Jeg vet ikke hva response etc bør være heller!?
fm <- fit(model)
summary(fm)
