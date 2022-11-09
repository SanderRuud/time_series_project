# Play with the data (start working on the project)
set.seed(1234)
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

# Hyperparameters and priors.
nstates <- 30 # Fixed number of states. 
pi <- rep(1, nstates)/nstates # Uniform prior state probabilities. 

# Try to fit a HMM. 
cell_train_df <- as.data.frame(cell_train)
#formel <- as.formula(paste0(paste(paste0("V", 1:71), collapse = "+"), "~ 1")) # Not sure if it is correct to simply add all the 71 time series?
# I think this is false, but not sure how to run one HMM on 71 different time series at once?!

# Quick way to make the list of formulas for the multivariate HMM. 
formel <- list()
families <- list()
for (i in 1:71){ # Only tried with 10 of the neurons for now, because of long fitting times. 
  # Running with all 71 neurons on Markov in order to see what happens. 
  formel[[i]] <- as.formula(paste0(paste0("V",i), "~1"))
  families[[i]] <- poisson(link = "log")
}

model <- depmix(formel, family = families, data = cell_train_df, 
                nstates = nstates, instart = pi) #  try to use pi instead for the initial state probabilities here. 
fm <- fit(model)
summary(fm)

#save(fm, file = "fitted_HMM.RData") # Save the fitted model such that we can load it locally.
load("fitted_HMM.RData", verbose = T)

# Now, what can we do with the states?
# How has Ben plotted them in the lecture on project description?
inferred_states <- fm@posterior[,1]
length(inferred_states)
table(inferred_states)

# Play with the transition matrix. 
transition_matrix <- matrix(getpars(fm)[(nstates+1):(nstates^2+nstates)],
                            byrow=TRUE,nrow=nstates,ncol=nstates)

#save(transition_matrix, file = "transition_matrix_full.RData") # Save the transition matrix such that we can load it locally.
load("transition_matrix_full.RData", verbose = T)

# Plot the transition matrix. 
transition_matrix %>% heatmap(Colv = NA)
transition_matrix %>% heatmap(Colv = NA, Rowv = NA) # Heatmap without the rearrangement because of the dendrogram.

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

# Try with a npn-linear dim. red. method, since PCA does not capture a lot of the variability in the first components. 
library(Rtsne)
tsne_out <- Rtsne(transition_matrix, pca = F, perplexity = as.integer(nstates/3)-1, theta = 0.8)
data.frame(cbind(tsne_out$Y,"state"=1:15)) %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  geom_text(aes(label = state), nudge_x = 1.5) +
  ggtitle(paste0("t-SNE on transition matrix")) +
  theme_minimal() # Not really a circle tbh. 


# Perhaps try to group the angular data in each state together, to see if each state "takes up" one specific direction.
# Then we could plot the recorded angles and see if they "coincide" with the states. 
hcl <- hclust(dist(transition_matrix), method = "ward.D")
cut <- cutree(hcl, k = 4) # We cut it into 4 groups. 
plot(hcl) # Produce the same dendrogram as in the heatmap above. 
rect.hclust(hcl, k = 4, border = 2:6)

# We group the inferred states according to the cut above. 
grouped_inferred_states <- cut[inferred_states]

# Group the angle data according to the states and then eventually by the cuts.
# Not sure how to do this right now?! HEEEEELP
