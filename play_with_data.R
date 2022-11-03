# Play with the data (start working on the project)

library(R.matlab)
library(ggplot2)
library(dplyr)
setwd("/home/ajo/gitRepos/time_series_project")
data <- readMat("Mouse28-140313_BS0150_HMMready.mat")
dim(data)
length(data)
str(data)

angdata <- data$resampledAwakeHeadAngleData
dim(angdata)
str(angdata)
sum(is.na(angdata)) # There are this many NaNs in the data?!
# This is the correct head direction angle data of the rats. This can be used for checking the performance of our model after fitting. 


celldata <- data$celldata
dim(celldata)
str(celldata)
# This is the celldata we are interested in using to infer the hidden states (which should represent the facing directions of the rats)

# Make train and test data for our model. Done "sequentially" in time to avoid data leakage.
N <- ncol(celldata)
train_indices <- as.integer(0.8*N) # The first 80% of observations are for training. 
cell_train <- t(celldata[,1:train_indices])
cell_test <- t(celldata[,(train_indices+1):N])
ang_train <- angdata[1:train_indices]
ang_test <- angdata[(train_indices+1):N]

# Fitting a HMM can e.g. be done with https://cran.r-project.org/web/packages/depmixS4/index.html
library(depmixS4)
?depmix

# Hyperparameters and priors.
nstates <- 10 # Fixed number of states. 
pi <- rep(1, nstates)/nstates # Uniform prior state probabilities. 

# Try to fit a HMM. 
cell_train_df <- as.data.frame(cell_train)
#formel <- as.formula(paste0(paste(paste0("V", 1:71), collapse = "+"), "~ 1")) # Not sure if it is correct to simply add all the 71 time series?
# I think this is false, but not sure how to run one HMM on 71 different time series at once?!

# Quick way to make the list of formulas for the multivariate HMM. 
formel <- list()
families <- list()
for (i in 1:10){ # Only tried with 10 of the neurons for now, because of long fitting times. 
  formel[[i]] <- as.formula(paste0(paste0("V",i), "~1"))
  families[[i]] <- poisson(link = "log")
}

model <- depmix(formel, family = families, data = cell_train_df, 
                nstates = nstates, instart = pi) #  try to use pi instead for the initial state probabilities here. 
fm <- fit(model)
summary(fm)

# Now, what can we do with the states?
# How has Ben plotted them in the lecture on project description?
inferred_states <- fm@posterior[,1]
length(inferred_states)
table(inferred_states)

# Play with the transition matrix. 
transition_matrix <- matrix(getpars(fm)[(nstates+1):(nstates^2+nstates)],
                            byrow=TRUE,nrow=nstates,ncol=nstates)

# Try to do PCA on transition matrix. 
pca <- princomp(transition_matrix)
summary(pca)
scores <- pca$scores
labels <- 1:nstates
df <- data.frame(cbind(scores[,1], scores[,2], labels))
df$labels<- as.factor(df$labels)
colnames(df) <- c("Z1", "Z2", "St")
plt <- tibble(df) %>% 
  ggplot(aes(x = Z1, y = Z2)) +
  geom_point() +
  geom_text(aes(label = St), nudge_x = 0.015) +
  ggtitle(paste0("PCA on transition matrix (first two pc's)")) +
  theme_minimal()
print(plt) # Not at all a circle :))


# Perhaps try to group the angular data in each state together, to see if each state "takes up" one specific direction.
# Then we could plot the recorded angles and see if they "coincide" with the states. 
