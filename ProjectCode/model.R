rm(list=ls())
set.seed(1234)
library(R.matlab)
library(ggplot2)
library(dplyr)
library(Rtsne)

data <- readMat("../Mouse28-140313_BS0150_HMMready.mat")
celldata <- data$celldata
celldata_df <- as.data.frame(t(celldata))
angdata <- data$resampledAwakeHeadAngleData
thresh <- rowSums(celldata)
THR <- 100 # I think this is a good threshold
celldata_active <- celldata_active[thresh > THR, ]
celldata_active_df <- as.data.frame(t(celldata_active))
nstates <- 35 # Fixed number of states. 
pies <- rep(1, nstates)/nstates # Uniform prior state probabilities. 


formel <- list()
families <- list()
n_neurons = 55 # Reduce this 
for (i in 1:n_neurons){ 
  formel[[i]] <- as.formula(paste0(paste0("V",i), "~1"))
  families[[i]] <- poisson()
}

# Part 1: 
model <- depmix(formel, family = families, data = celldata_df, 
                nstates = nstates, instart = pies) #  try to use pi instead for the initial state probabilities here. 
fm <- fit(model)

write.csv((fm@posterior), "../posterior.csv")
nstates = length(fm@posterior)-1
transition_matrix <- matrix(getpars(fm)[(nstates+1):(nstates^2+nstates)],
                            byrow=TRUE,nrow=nstates,ncol=nstates)
write.csv(transition_matrix, "../trans_matrix.csv")



# Part 2: Only considering active neurons

model2 <- depmix(formel, family = families, data = celldata_active_df, 
                nstates = nstates, instart = pies) #  try to use pi instead for the initial state probabilities here. 
fm2 <- fit(model)

write.csv((fm@posterior), "../posterior2.csv")
nstates = length(fm@posterior)-1
transition_matrix2 <- matrix(getpars(fm)[(nstates+1):(nstates^2+nstates)],
                            byrow=TRUE,nrow=nstates,ncol=nstates)
write.csv(transition_matrix2, "../trans_matrix2.csv")


# Generating plots

inferred_states <- fm@posterior[,1]
length(inferred_states)
table(inferred_states)

# Play with the transition matrix. 

# Plot the transition matrix. 
pdf(file = "trans_matrix.pdf", width = 9, height = 5)
transition_matrix %>% heatmap(Rowv = NA)
dev.off()
pdf(file = "trans_matrix2.pdf", width = 9, height = 5)
transition_matrix %>% heatmap(Colv = NA, Rowv = NA) # Heatmap without the rearrangement because of the dendrogram.
dev.off()

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

pdf(file = "PCA.pdf", width = 9, height = 5)
print(plt) # Not at all a circle :))
dev.off()

# Try with a npn-linear dim. red. method, since PCA does not capture a lot of the variability in the first components. 
tsne_out <- Rtsne(transition_matrix, pca = F, perplexity = as.integer(30/3)-1, theta = 0.8)
t_sne <- data.frame(cbind(tsne_out$Y,"state"=1:15)) %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  geom_text(aes(label = state), nudge_x = 1.5) +
  ggtitle(paste0("t-SNE on transition matrix")) +
  theme_minimal() # Not really a circle tbh. 

pdf(file = "t-SNE.pdf", width = 9, height = 5)
print(t_sne)
dev.off()

# Perhaps try to group the angular data in each state together, to see if each state "takes up" one specific direction.
# Then we could plot the recorded angles and see if they "coincide" with the states. 
hcl <- hclust(dist(transition_matrix), method = "ward.D")
cut <- cutree(hcl, k = 4) # We cut it into 4 groups. 

pdf(file = "hclust.pdf", width = 9, height = 5)
plot(hcl) # Produce the same dendrogram as in the heatmap above. 
rect.hclust(hcl, k = 4, border = 2:6)
dev.off()

# We group the inferred states according to the cut above. 
grouped_inferred_states <- cut[inferred_states]
