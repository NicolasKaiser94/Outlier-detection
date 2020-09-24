if (!require("sampling")) install.packages("sampling")
if (!require("MASS")) install.packages("MASS")
if (!require("cluster")) install.packages("cluster", "factoextra")
if (!require("PMCMRplus")) install.packages("PMCMRplus")
if (!require("kohonen")) install.packages("kohonen")
if (!require("OutlierDetection"))install.packages("OutlierDetection")

library(sampling)
library(MASS)
library(cluster)
library(factoextra)
library(PMCMRplus)
library(kohonen)
library(OutlierDetection)
library(ggplot2)
library(gridExtra)



# Defining the size of our simple random sample without replacement
n <- 1000

# Converting our Incomevariable into a dataframe
INC <- as.data.frame(INC)

# Drawing a Simple random sample from the synthetic AMELIA Dataset. This serves as a new 

SI_SAMPLE <- srswor(n, nrow(INC))
Population <- INC[SI_SAMPLE>0, ] # Perfect simple random sample without missing data

# Comparing population and sample
summary(INC)
summary(Population)

# Generating a density curve

plot(density(Population))

## Comparing outliers with simple boxplots
boxplot(Population)

# Using the various for outlier detection in detecting possible outliers

## Z score method
z_scores <- scale(Population)
Zscore_outliers <- Population[which(scale(Population) > 3)]

### Generalized Extreme Studentized Deviate 

gesd_res <- gesdTest(Population, maxr = 500)

outliers.GESD <- Population[which(gesd_res$p.value < 0.05)] 

which(gesd_res$p.value < 0.05)


## Inter Quartile range method

Q1<-quantile(Population,0.25)
Q3<-quantile(Population,0.75)
IQR<-(Q3-Q1)

IQR.outliers <- function(x) {
  Q3<-quantile(x,0.75, na.rm = TRUE)
  Q1<-quantile(x,0.25, na.rm = TRUE)
  IQR<-(Q3-Q1)
  lower.bound<- (Q1-(1.5*IQR))
  upper.bound<- (Q3+(1.5*IQR))
  c(x[x <lower.bound],x[x>upper.bound])

  }
IQR.outliers(Population)


#Elbow Method for finding the optimal number of clusters
set.seed(2024)

# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- scale(Population)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")



## k Nearest Neighbour

K.outlier <-nn(as.matrix(Population), k = 3, cutoff = 0.95, Method = "euclidean",
   rnames = FALSE, boottimes = 1000)



## K Medoids Method

set.seed(2024)

kmedoids_res <- clara(Population, 3)

plot(kmedoids_res)

kmedoids_res$medoids

kmed.centers <- kmedoids_res$medoids[kmedoids_res$clustering, ]

kmed.distance <- sqrt(((Population - kmed.centers)^2))

kmedoids_out <- order(kmed.distance, decreasing=T)

kmed.outliers <- head(kmedoids_out, 20) # represent the first 20 large distance seen to be outliers after plotting. 
kmedoids_outliers <- Population[kmed.outliers] # Data points seen as outliers

plot(Population, pch=16, col="blue")
points(kmed.outliers, Population[kmed.outliers], pch=16, col="red") 


### Self Organising Maps (unsupervised)
set.seed(2024)

X_som <- data.frame(Population)
X_som <- scale(X_som)

som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")

som_model <- som(X_som, 
                 grid=som_grid, 
                 rlen=1000, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE, radius= 1 )

summary(som_model)


## Visualisation

plot(som_model, type="changes") 

plot(som_model, type="count")              # Node counts plot

plot(som_model, type="dist.neighbours")    # Neighbour distance plot

plot(som_model, type="codes", palette.name = rainbow)   # Codes / Weight vectors

plot(som_model, type="mapping") # Visualisation of th heat map


som_distance <- order(som_model$distances, decreasing = T)

som_out <-head(som_distance, 100)

som_outliers <- Population[som_out]



# 6- We then simulate two separate  univariate normal data A and B using the mean and covariance from 
# the "clean" and "contaminated" data.

# 7. Next we randomly sample some data points from B with some proportions like (0, 0.01, 0.025, 0.05)
# and using this equation => (1-p)A + pB where p is the proportion .

# 8.	We then add those data point to A. (In the templ paper they swapped a cell in data A with a 
# cell from data B). By doing this we know the outliers we are adding

# 9. And then we apply the methods again to the data that will result from (7) and see how many of the 
# ouliers our methods can "rightly" detect.


#############
#
# Monte-Carlo Simulation for simulation study. Create a function
# Before we create our function, we check if all the single steps work properly.

SI_SAMPLE <- srswor(500000, nrow(INC)) ### Sampling from the synthetic INCome distribution.
Population <- INC[SI_SAMPLE>0, ] ### It serves as our population
SI_SAMPLE_2 <- srswor(250, nrow(as.data.frame(Population))) ### Sampling from the population to obtain a sample
Sample <- as.data.frame(Population)[SI_SAMPLE_2>0, ]

z_scores <- scale(Sample) # Applying the first Method to evaluate: The z-score method
Zscore_outliers <- Sample[which(scale(Sample) > 3)] # Those are the outliers detected by the method

gesd_res <- gesdTest(Sample, maxr = 20) # Applying the GESD
outliers.GESD <- Sample[which(gesd_res$p.value < 0.05)] 

IQR_outliers <- IQR.outliers(Sample) # Applying the IQR-method

#Finding the optimal number of k clusters using kmeans
k.max <- 15    
data <- scale(Sample)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})



## K nearest neighbour
knearest.out <- nn(as.matrix(Sample), k= 3,cutoff = 0.95, Method = "euclidean",
                   rnames = FALSE, boottimes = 1000)

knn_out <- knearest.out$`Outlier Observations`


## K-medoids

kmedoids_res <- clara(Population, 3)


kmed.centers <- kmedoids_res$medoids[kmedoids_res$clustering, ]

kmed.distance <- sqrt(((Population - kmed.centers)^2))

kmedoids_out <- order(kmed.distance, decreasing=T)

kmed.outliers <- head(kmedoids_out, 20) # represent the first 20 large distance seen to be outliers after plotting.

kmedoids_outliers <- Population[kmed.outliers] # Data points seen as outliers


## Self Organising Maps

X_som <- data.frame(Sample)
X_som <- scale(X_som)

som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")

som_model <- som(X_som, 
                 grid=som_grid, 
                 rlen=1000, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE, radius= 1 )

som_distance <- order(som_model$distances, decreasing = T)

som_out <-head(som_distance, 100)

som_outliers <- Population[som_out]




# We divide the dataset into a "clean" part (no outliers) and a "contaminated" part (all outliers)
# So far, "true" outliers are all the outliers which were found by at least one method. (Methodological questionable?!)
#contaminated <- Sample[which(Sample %in% Zscore_outliers & Sample %in% outliers.GESD|
                             #  Sample %in% Zscore_outliers & Sample %in% IQR_outliers|
                              # Sample %in% outliers.GESD & Sample %in% IQR_outliers)] # Dividing the datatset into
contaminated <- Sample[which(Sample %in% Zscore_outliers|Sample %in% outliers.GESD|Sample %in% IQR_outliers|
                               Sample %in% knn_out|Sample %in% kmedoids_outliers |Sample %in% som_outliers)]
clean <- Sample[!(Sample %in% contaminated)] # a "clean" and a "contaminated" dataset

clean_mean <-  mean(clean) # Saving the mean and standard deviation from the "clean" and "contaminated"
contaminated_mean <- mean(contaminated) # Dataset. This wil be needed for the generation of
sdev_clean <- sd(clean)                 # The two normal distributed data sets
sdev_cont <- sd(contaminated)

new_data_a <- rnorm(n = 250, mean = clean_mean, sd = sdev_clean) # Creating distribution based on the "clean" data
new_data_b <- rnorm(n = 1000, mean = contaminated_mean, sd = sdev_cont) # Creating a distribution based on the "contaminated" data

outlier_n <- 20 # Number of outliers

sampled_outliers_process <- srswor(n = outlier_n, nrow(as.data.frame(new_data_b))) # Sampling outliers from distribution B
sampled_outliers <- as.data.frame(new_data_b)[sampled_outliers_process>0,]

new_data_a_b <- data.frame(rbind(as.matrix(new_data_a), as.matrix(sampled_outliers))) # Adding the sampled outliers to the new generated data
# The rbind command appends the outliers to the end of the data frame. Therefore, we know
# the "true" and "false outliers.
plot(density(new_data_a_b[,1]))
false <- new_data_a_b[c(nrow(new_data_a_b) - outlier_n : nrow(new_data_a_b)),] 
true <- new_data_a_b[-c(nrow(new_data_a_b) - outlier_n : nrow(new_data_a_b)),]


## Now we apply the methods again. We differentiate between correctly detected outliers
# and falsely detected outliers (false alarm rate)
# Z-score method

z_scores <- scale(new_data_a_b)
Zscore_outliers <- new_data_a_b[which(scale(new_data_a_b) > 3),]

gesd_res <- gesdTest(as.matrix(new_data_a_b), maxr = 20) # Applying the GESD
outliers.GESD <- new_data_a_b[which(gesd_res$p.value < 0.05),] 

IQR.out <- IQR.outliers(as.matrix(new_data_a_b))

k.max <- 15
data <- scale(new_data_a_b[,1])
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss}) 
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") # From the graph k = 3


knearest.out <- nn(as.matrix(new_data_a_b), k= 3, boottimes = 1000) #Applying k-nearest neighbor

knn_out <- knearest.out$`Outlier Observations`


kmedoids_res <- clara(new_data_a_b, 3) #Applying k medoids


kmed.centers <- kmedoids_res$medoids[kmedoids_res$clustering, ]

kmed.distance <- sqrt(((new_data_a_b - kmed.centers)^2))

kmedoids_out <- order(kmed.distance, decreasing=T)

kmed.outliers <- head(kmedoids_out, 20)

kmedoids_outliers <- new_data_a_b[kmed.outliers,]


## Sef organising maps

X_som <- data.frame(new_data_a_b)
X_som <- scale(X_som)

som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")

som_model <- som(X_som, 
                 grid=som_grid, 
                 rlen=1000, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE, radius= 1 )

som_distance <- order(som_model$distances, decreasing = T)

som_out <-head(som_distance, 50)

som_outliers <- new_data_a_b[som_out,]


length(intersect(Zscore_outliers, true)) / length(true) * 100 # shows outliers common in both
length(intersect(Zscore_outliers, false)) / length(false) * 100

length(intersect(IQR.out, true)) / length(true) * 100 # shows outliers common in both
length(intersect(IQR.out, false)) / length(false) * 100

length(intersect(outliers.GESD, true)) / length(true) * 100 # shows outliers common in both
length(intersect(outliers.GESD, false)) / length(false) * 100

length(intersect(knn_out, true)) / length(false) * 100
length(intersect(knn_out, false)) / length(false) * 100

length(intersect(kmedoids_outliers, true)) / length(false) * 100
length(intersect(kmedoids_outliers, false)) / length(false) * 100


length(intersect(som_outliers, true)) / length(false) * 100
length(intersect(som_outliers, false)) / length(false) * 100


### Now we implement this as a function, which allows us to vary the relevant variables as well as
### to perform our Monte-Carlo-Simulation
outlier_simulation <- function(N = N, n = n, outlier_n = outlier_n,
                               dist = dist, Z_score_threshold = 
                               Z_score_threshold, GESD_p_value = GESD_p_value, k_nearest_k,
                               new_data_a_n = new_data_a_n, new_data_b_n = new_data_a_n,
                               maxr = maxr, k_medoids_outliers = k_medoids_outliers, som_maps_outliers = som_maps_outliers) 
{
  SI_SAMPLE <- srswor(N, nrow(INC)) ### Sampling from the synthetic INCome distribution.
  Population <- INC[SI_SAMPLE>0, ] ### It serves as our population
  SI_SAMPLE_2 <- srswor(n, nrow(as.data.frame(Population))) ### Sampling from the population to obtain a sample
  Sample <- as.data.frame(Population)[SI_SAMPLE_2>0, ]
  
  z_scores <- scale(Sample) # Applying the first Method to evaluate: The z-score method
  Zscore_outliers <- Sample[which(scale(Sample) > Z_score_threshold)] # Those are the outliers detected by the method
  
  gesd_res <- gesdTest(Sample, maxr = maxr) # Applying the GESD
  outliers.GESD <- Sample[which(gesd_res$p.value < GESD_p_value)] 
  
  IQR_outliers <- IQR.outliers(Sample) # Applying the IQR-method

  
  knearest.out <- nn(as.matrix(Sample),cutoff = 0.95, Method = "euclidean", # Aoolying K-Nearest Neighbor
                     rnames = FALSE, k = k_nearest_k, boottimes = 1000)
  
  knn_out <- knearest.out$`Outlier Observations`
  
  
  kmedoids_res <- clara(Sample, 3)   ### Applying K Medoids
  
  
  kmed.centers <- kmedoids_res$medoids[kmedoids_res$clustering, ]
  
  kmed.distance <- sqrt(((Sample - kmed.centers)^2))
  
  kmedoids_out <- order(kmed.distance, decreasing=T)
  
  kmed.outliers <- head(kmedoids_out, k_medoids_outliers)
  
  kmedoids_outliers <- Sample[kmed.outliers]
  
  
  X_som <- data.frame(Sample)         #### Applying Self Organising Maps
  X_som <- scale(X_som)
  
  som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")
  
  som_model <- som(X_som, 
                   grid=som_grid, 
                   rlen=1000, 
                   alpha=c(0.05,0.01), 
                   keep.data = TRUE, radius= 1 )
  
  som_distance <- order(som_model$distances, decreasing = T)
  
  som_out <-head(som_distance, som_maps_outliers)
  
  som_outliers <- Sample[som_out]
  
  
  #contaminated <- Sample[which(Sample %in% Zscore_outliers & Sample %in% outliers.GESD|
  #                               Sample %in% Zscore_outliers & Sample %in% IQR_outliers|
  #                               Sample %in% outliers.GESD & Sample %in% IQR_outliers)] # Dividing the datatset into
  contaminated <- Sample[which(Sample %in% Zscore_outliers|Sample %in% outliers.GESD|Sample %in% IQR_outliers|
                                 Sample %in% knn_out|Sample %in% kmedoids_outliers|Sample %in% som_outliers)]
  
  clean <- Sample[!(Sample %in% contaminated)] # a "clean" and a "contaminated" dataset
  
  clean_mean <-  mean(clean) # Saving the mean and standard deviation from the "clean" and "contaminated"
  contaminated_mean <- mean(contaminated) # Dataset. This wil be needed for the generation of
  sdev_clean <- sd(clean)                 # The two normal distributed data sets
  sdev_cont <- sd(contaminated)
 
  new_data_a <- rnorm(n = new_data_a_n, mean = clean_mean, sd = sdev_clean) # Creating distribution based on the "clean" data
  new_data_b <- rnorm(n = new_data_b_n, mean = contaminated_mean, sd = sdev_cont) # Creating a distribution based on the "contaminated" data

  sampled_outliers_process <- srswor(n = outlier_n, nrow(as.data.frame(new_data_b))) # Sampling outliers from distribution B
  sampled_outliers <- as.data.frame(new_data_b)[sampled_outliers_process>0,]
  
  new_data_a_b <- data.frame(rbind(as.matrix(new_data_a), as.matrix(sampled_outliers))) # Adding the sampled outliers to the new generated data
  # The rbind command appends the outliers to the end of the data frame. Therefore, we know
  # the "true" and "false outliers.
  
  false <- new_data_a_b[c(nrow(new_data_a_b) - outlier_n : nrow(new_data_a_b)),] 
  true <- new_data_a_b[-c(nrow(new_data_a_b) - outlier_n : nrow(new_data_a_b)),]
  
  z_scores <- scale(new_data_a_b) # Applying the Zscore
  Zscore_outliers <- new_data_a_b[which(scale(new_data_a_b) > 3),]
  
  gesd_res <- gesdTest(as.matrix(new_data_a_b), maxr = maxr) # Applying the GESD
  outliers.GESD <- new_data_a_b[which(gesd_res$p.value < 0.05),] 
  new_data_a_b <- as.matrix(new_data_a_b)
  IQR.out <- IQR.outliers(new_data_a_b) # Applying the IQR
  
  knearest.out <- nn(as.matrix(new_data_a_b), k= k_nearest_k, cutoff = 0.95, Method = "euclidean",
                     rnames = FALSE, boottimes = 1000)             #Applying the K-nearest neighbor
  knn_out <- knearest.out$`Outlier Observations` 
 
  kmedoids_res <- clara(new_data_a_b[,1], 3)
  
  
  kmed.centers <- kmedoids_res$medoids[kmedoids_res$clustering, ] # Applying K Medoids
  
  kmed.distance <- sqrt(((new_data_a_b[,1] - kmed.centers)^2))
  
  kmedoids_out <- order(kmed.distance, decreasing=T)
  
  kmed.outliers <- head(kmedoids_out, k_medoids_outliers)
  
  kmedoids_outliers <- new_data_a_b[,1][kmed.outliers]
  
  
  X_som <- data.frame(new_data_a_b[,1])          #### Applying Sef Organising Maps
  X_som <- scale(X_som)
  
  som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")
  
  som_model <- som(X_som, 
                   grid=som_grid, 
                   rlen=1000, 
                   alpha=c(0.05,0.01), 
                   keep.data = TRUE, radius= 1 )
  
  som_distance <- order(som_model$distances, decreasing = T)
  
  som_out <-head(som_distance, som_maps_outliers)
  
  som_outliers <- new_data_a_b[,1][som_out]
  
  
  True_zscore <- length(intersect(Zscore_outliers, true)) / length(true) * 100 # shows outliers common in both
  False_zscore <- length(intersect(Zscore_outliers, false)) / length(false) * 100
  
  True_IQR <- length(intersect(IQR.out, true)) / length(true) * 100 # shows outliers common in both
  False_IQR <- length(intersect(IQR.out, false)) / length(false) * 100
  
  True_GESD <- length(intersect(outliers.GESD, true)) / length(true) * 100 # shows outliers common in both
  False_GESD <- length(intersect(outliers.GESD, false)) / length(false) * 100
  
  True_k_neighbor <- length(intersect(knn_out, true)) / length(true) * 100
  False_k_neighbor <- length(intersect(knn_out, false)) / length(false) * 100

  True_k_medoids <- length(intersect(kmedoids_outliers, true)) / length(true) * 100
  False_k_medoids <- length(intersect(kmedoids_outliers, false)) / length(false) * 100
  
  True_som <- length(intersect(som_outliers, true)) / length(true) * 100
  False_som <- length(intersect(som_outliers, false)) / length(false) * 100
    
  mylist <- list(True_zscore = True_zscore, False_zscore = False_zscore, True_IQR = True_IQR,
                 False_IQR = False_IQR, True_GESD = True_GESD, False_GESD = False_GESD,
                 True_k_neighbor = True_k_neighbor, False_k_neighbor = False_k_neighbor,
                 True_k_medoids = True_k_medoids, False_k_medoids = False_k_medoids,
                 True_som = True_som, False_som = False_som)
  
  return(mylist)
}
outlier_n <- 100
outlier_simulation(N = 500000, n = 10000, outlier_n = outlier_n, Z_score_threshold = 3, GESD_p_value = 0.05,
                   k_nearest_k = 3,
                   new_data_a_n = n, new_data_b_n = 10000, maxr = 5, 
                   k_medoids_outliers = outlier_n, 
                   som_maps_outliers = outlier_n)
                   

## Now we use a for-loop to apply a Monte-carlo simulation
outlier_n <- 10
set.seed(1994)
result <- list() # Results of each simulation run will be stored here
for(i in 1:1000){ # Number of Simulation runs
  resultX <- outlier_simulation(N = 500000, # Number of elements in synthetic population
                                     n = 1000,   # Sample size (SRS)
                                     outlier_n =  outlier, # Number of outliers sampled from contaminated data set
                                     k_nearest_k = 3 #
                                     Z_score_threshold = 3, # Threshold for Z-score method
                                     GESD_p_value = 0.05, # p-value threshold for GESD                                     
                                     new_data_a_n = n-outlier_n, # Size of "clean" data set
                                     new_data_b_n = 2000, # Size of "contaminated" data set
                                     maxr = outlier_n, # Maximum number of outliers GESD 
                                    k_medoids_outliers = outlier_n, # Number of data points considered outliers
                                    som_maps_outliers= outlier_n) # Number of data points considered as outliers
                                      
                                      
  result[[length(result)+1]] <- resultX # Each simulation run is stored in a list
  print(i)
}

# Converting the results into a readable form

Res_Norm_1 <- c((True_Zscore_result <- mean(unlist(lapply(result, "[", 1)))),
(False_Zscore_result <- mean(unlist(lapply(result, "[", 2)))),
(True_IQR_result <- mean(unlist(lapply(result, "[", 3)))),
(False_IQR_result <- mean(unlist(lapply(result, "[", 4)))),
(True_GESD_result <- mean(unlist(lapply(result, "[", 5)))),
(False_GESD_result <- mean(unlist(lapply(result, "[", 6)))),
(True_k_neighbor_result <- mean(unlist(lapply(result, "[", 7)))),
(False_k_neighbor_result <- mean(unlist(lapply(result, "[", 8)))),
(True_k_medoids_result <- mean(unlist(lapply(result, "[", 9)))),
(False_k_medoids_result <- mean(unlist(lapply(result, "[", 10)))),
(True_som_result <- mean(unlist(lapply(result, "[", 11)))),
(False_som_result <- mean(unlist(lapply(result, "[", 12)))))


  


Res_Norm_1  # (1000, 0.01) number of outlier = 10
Res_Norm_2 # (1000, 0.025) number of outlier = 25
Res_Norm_3 # (1000, 0.05)  number of outlier = 50

Res_Norm_4  # (10000, 0.01) number of outliers = 100
Res_Norm_5 # (10000, 0.025) number of outliers = 250
Res_Norm_6 # (10000, 0.05)  number of outliers = 500

Res_Norm_4  # (100000, 0.01) number of outliers = 1000
Res_Norm_5 # (100000, 0.025) number of outliers = 2500
Res_Norm_6 # (100000, 0.05)  number of outliers = 5000











setwd("C:/Users/User/Desktop/M.sc. Applied Statistics/3. Semester/Research Project")
simulation_results <- read.csv("Simulation_results_Norm.csv", header = TRUE)
simulation_results <- read.csv("Simulation_results_Exponential.csv", header = TRUE)

# Converting Dataset size and Percentage of outliers to factor.

simulation_results$Percentage<- as.factor(simulation_results$Percentage)
simulation_results$Size <- as.factor(simulation_results$Size)

# Converting factors to numeric values
simulation_results$Detection.Rate <- gsub(",", "", simulation_results$Detection.Rate)
simulation_results$Detection.Rate <- as.numeric(simulation_results$Detection.Rate)
simulation_results$Detection.Rate <- simulation_results$Detection.Rate / 100

simulation_results$False <- gsub(",", "", simulation_results$False)
simulation_results$False <- as.numeric(simulation_results$False)
simulation_results$False <- simulation_results$False / 100

str(simulation_results)

par(mfrow=c(1,3))

figure_1 <- ggplot(data = simulation_results[which(simulation_results$Size == 1000),],
                   aes(x = Percentage, y = Detection.Rate, 
                       col = Method,
                       group = Method)) + 
  geom_point() +
  geom_line() +
  ylim(c(0, 75)) +
  xlab("Amount of contamination in %") +
  ylab("Detection rate in %") +
  theme_classic()
figure_1
figure_2 <- ggplot(data = simulation_results[which(simulation_results$Size == 10000),],
                   aes(x = Percentage, y = Detection.Rate, 
                       col = Method,
                       group = Method)) + 
  geom_point() +
  geom_line() +
  ylim(c(0, 75)) +
  xlab("Amount of contamination in %") +
  ylab("Detection rate in %") +
  theme_classic()
figure_2

figure_3 <- ggplot(data = simulation_results[which(simulation_results$Size == 20000),],
                   aes(x = Percentage, y = Detection.Rate, 
                       col = Method,
                       group = Method)) + 
  geom_point() +
  geom_line() +
  ylim(c(0, 75)) +
  xlab("Amount of contamination in %") +
  ylab("Detection rate in %") +
  theme_classic()
figure_3

grid.arrange(figure_1, figure_2, figure_3, ncol=3)
figure_4 <- ggplot(data = simulation_results[which(simulation_results$Size == 1000),],
                   aes(x = Percentage, y = False, 
                       col = Method,
                       group = Method)) + 
  geom_point() +
  geom_line() +
  ylim(c(0, 5)) +
  xlab("Amount of contamination in %") +
  ylab("False alarm rate in %") +
  theme_classic()
figure_4
figure_5 <- ggplot(data = simulation_results[which(simulation_results$Size == 10000),],
                   aes(x = Percentage, y = False, 
                       col = Method,
                       group = Method)) + 
  geom_point() +
  geom_line() +
  ylim(c(0, 5)) +
  xlab("Amount of contamination in %") +
  ylab("False alarm rate in %") +
  theme_classic()
figure_5
figure_6 <- ggplot(data = simulation_results[which(simulation_results$Size == 20000),],
                   aes(x = Percentage, y = False, 
                       col = Method,
                       group = Method)) + 
  geom_point() +
  geom_line() +
  ylim(c(0, 5)) +
  xlab("Amount of contamination in %") +
  ylab("False alarm rate in %") +
  theme_classic()
figure_6
grid.arrange(figure_1, figure_2, figure_3, figure_4, figure_5, figure_6, ncol=3)
