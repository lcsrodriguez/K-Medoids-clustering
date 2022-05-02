## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
graphics.off()


## ---- results='hide'----------------------------------------------------------------------------------------------------------------------------------
df <- iris # Dataset d'origine
df_init <- df[1:4] # Dataset sans la dernière colonne (on supprime la dernière colonne)
X <- df_init
k <- 3
algo_k_medoids <- function(X, k){
  # Récupération des dimensions 
  n <- nrow(X)
  p <- ncol(X)
  # Etape 1.1  =================================================
  d <- data.frame(matrix(data = 0, nrow = n, ncol = n))
  colnames(d) <- 1:n
  
  # Calcul de d
  for (i in 1:n) {
    for (j in 1:n) {
      s <- 0
      for (a in 1:p) {
        s <- s + (X[i, a] - X[j, a])^2
      }
      d[i, j] <- sqrt(s)
    }
  }
  
  # Etape 1.2  =================================================
  v <- data.frame(matrix(data = 0, nrow = n, ncol = 2))
  colnames(v) <- c("objet", "v")
  v$objet <- 1:n
  
  # Calcul de v_j pour chaque j
  for (j in 1:n) {
    
    v_j <- 0
    for (i in 1:n) {
      
      denom <- 0
      for (l in 1:n) {
        denom <- denom + d[i, l]
      }
      v_j <- v_j + d[i, j]/denom
    }
    v[j, 2] <- v_j
  }
  
  # Etape 1.3  =================================================
  v_sorted <- v[order(v[,ncol(v)], decreasing = FALSE),]
  # Sélection des k objets en tant que médoides initiaux
  medoids <- v_sorted[1:k,]$objet
  
  # Etape 1.4  =================================================
  # Création du cluster initial
  clusters <- data.frame(matrix(data = 0, nrow = n, ncol = 1))
  clusters_reals <- data.frame(matrix(data = 0, nrow = n, ncol = 1))
  colnames(clusters) <- c("cluster")
  colnames(clusters_reals) <- c("c")
  for (i in 1:n) { # Pour chaque objet i dans 1:n, on lui associe la classe dont le médoide est le plus proche (minimisation distance objet <-> médoide)
    # Initialisation de la matrice des clusters
    clusters[i, 1] <- -1 # Initialisation, au départ, aucun objet n'appartient à une classe
    
    
    d_min <- d[medoids[1], i]
    # Détermination de la classe pour chaque objet
    for (j in 1:k) {
      if (d[i, medoids[j]] <= d_min ) {
        clusters[i, 1] <- medoids[j]
        clusters_reals[i, 1] <- j
        d_min <- d[i, medoids[j]]
      }
    }
  }
  
  # Etape 1.5  =================================================
  # Calcul de la somme des distances entre tous les objets et leurs médoïdes respectifs
  # Remarque : Pour obtenir le médoîde d'un objet i, on prend clusters[, 2][i]
  sum_distances <- 0
  for (i in 1:n) {
    sum_distances <- sum_distances + d[i, clusters$cluster[i]]
  }
  
  old_sum_distances <- sum_distances
  new_sum_distances <- old_sum_distances * 2 # Initialisation de cette quantité au double de la première pour garantir le premier passage dans la boucle
  epsilon <- 1e-1
  compteur_iter <- 0
  ITER_MAX <- 26
  while (abs(old_sum_distances - new_sum_distances) > epsilon & compteur_iter <= ITER_MAX ) {
    #print("############ Début itération")
    old_sum_distances <- new_sum_distances
    
    # Etape 2  ===================================================
    # Pour chaque cluster
    for (l in 1:k) {
      
      # Récupération des objets de chaque cluster
      objets_cluster <- as.numeric(rownames(clusters)[clusters$cluster == medoids[l]])
      
      # Récupération médoïde actuel du cluster
      current_medoid <- medoids[l]
      #print("Médoïde actuel :")
      #print(current_medoid)
    
      #print("Objets clusters : ")
      #print(objets_cluster)
      #print("==================")
      
      # Calcul des distances totales (total distance to others objects)
      N <- length(objets_cluster)
      total_distances <- data.frame(matrix(data = 0, ncol = 1, nrow = N))
      rownames(total_distances) <- objets_cluster
      colnames(total_distances) <- c("distance")
      
      for (q in objets_cluster) {
        tt_dist <- 0
        for (r in objets_cluster) {
          tt_dist <- tt_dist + d[q, r]
        }
        #print(tt_dist)
        total_distances$distance[rownames(total_distances) == q] <- tt_dist
      }
      
      # Détermination nouveau médoïde (objet minimisant le tableau des distances totales)
      new_medoid <- which.min(total_distances$distance)
      #print("Nouveau médoïde")
      #print(new_medoid)
      
      # Remplacement nouveau médoïde
      medoids[l] <- new_medoid
    }
    # Etape 3.1  =================================================
    # Attribuer chaque objet au médoïde le plus proche (étape semblable à l'étape 1.4)
    for (i in 1:n) { # Pour chaque objet i dans 1:n, on lui associe la classe dont le médoide est le plus proche (minimisation distance objet <-> médoide)
      d_min <- d[medoids[1], i]
      # Détermination de la classe pour chaque objet
      for (j in 1:k) {
        if (d[i, medoids[j]] <= d_min) {
          clusters[i, 1] <- medoids[j]
          clusters_reals[i, 1] <- j
          d_min <- d[i, medoids[j]]
        }
      }
    }
    
    # Etape 1.5  ================================================= (calcul de new_sum_distances)
    # Calcul de la somme des distances entre tous les objets et leurs médoïdes respectifs
    
    # Remarque : Pour obtenir le médoîde d'un objet i, on prend clusters[, 2][i]
    sum_distances <- 0
    for (i in 1:n) {
      sum_distances <- sum_distances + d[i, clusters$cluster[i]]
    }
    new_sum_distances <- sum_distances
    
    #print("############ Fin itération (print old puis new)")
    #print(old_sum_distances)
    #print(new_sum_distances)
    compteur_iter <- compteur_iter + 1
  }
  print("Fin algorithme")
  print("Nb d'itérations : ")
  print(compteur_iter)
  return(list(clusters = clusters_reals$c, medoids = clusters))
}
#res_k_medoids_maison <- algo_k_medoids(X, k)$clusters


## -----------------------------------------------------------------------------------------------------------------------------------------------------
n <- 120 # Nombre d'objets
simulation_dataset_cluster_A <- function(){
  # Simulation du cluster A
  mu_A_x <- 0
  mu_A_y <- 0
  sigma_A <- 1.5
  cluster_A <- cbind.data.frame(x = rnorm(n%/%3, mean = mu_A_x, sd = sigma_A), y = rnorm((n%/%3), mean = mu_A_y, sd = sigma_A))
  return(cluster_A)
}

simulation_dataset_cluster_B <- function(){
  # Simulation du cluster B
  mu_B_x <- 6
  mu_B_y <- -1
  sigma_B <- 0.5
  cluster_B <- cbind.data.frame(x = rnorm(n%/%3, mean = mu_B_x, sd = sigma_B), y = rnorm((n%/%3), mean = mu_B_y, sd = sigma_B))
  return(cluster_B)
}

simulation_dataset_cluster_C <- function(){
  # Simulation du cluster C
  mu_C_x <- 6
  mu_C_y <- 2
  sigma_C <- 0.5
  sigma_C_L <- 2
  proportion <- 0.1
  a <- rnorm(proportion*(n%/%3), mean = mu_C_x, sd = sigma_C_L)
  cluster_C_L <- cbind.data.frame(x = rnorm(proportion*(n%/%3), mean = mu_C_x, sd = sigma_C_L), y = rnorm(proportion*(n%/%3), mean = mu_C_y, sd = sigma_C_L))
  cluster_C_2 <- cbind.data.frame(x = rnorm((1 - proportion)*(n%/%3), mean = mu_C_x, sd = sigma_C), y = rnorm((1 - proportion)*(n%/%3), mean = mu_C_y, sd = sigma_C))
  cluster_C <- rbind(cluster_C_L, cluster_C_2)
  return(cluster_C)
}

cluster_A <- simulation_dataset_cluster_A()
cluster_B <- simulation_dataset_cluster_B()
cluster_C <- simulation_dataset_cluster_C()

# Jointure afin de tester sur les algorithmes
cluster_test <- rbind(cluster_A, cluster_B, cluster_C)

simulation_dataset <- function(){
  cluster_A <- simulation_dataset_cluster_A()
  cluster_B <- simulation_dataset_cluster_B()
  cluster_C <- simulation_dataset_cluster_C()

  # Jointure afin de tester sur les algorithmes
  cluster_test <- rbind(cluster_A, cluster_B, cluster_C)
  return(cluster_test)
}

# Affichage graphique des clusters
library(ggplot2)
ggplot(cluster_A, aes(x = x, y = y), main = "d") + 
  geom_point(colour = "red") + 
  geom_point(data = cluster_B, colour = "green") + 
  geom_point(data = cluster_C, colour = "blue") + 
  ggtitle("Représentation graphique du dataset artificiellement généré") + 
  xlab("Abscisse (x)") + 
  ylab("Ordonnée (y)")


## ---- results = "hide"--------------------------------------------------------------------------------------------------------------------------------
library(fossil)
library(cluster)
library(mclust)
n_sim <- 20
k <- 3

rands <- data.frame(matrix(0, nrow = n_sim, ncol = 4))
colnames(rands) <- c("maison", "pam", "kmeans", "mclust")


for (sim in 1:n_sim) {
  # Simulation du dataset
  cluster_test <- scale(simulation_dataset())
  cluster_test_real <- rbind(data.frame(c = rep(1, (n%/%3))), data.frame(c = rep(2, (n%/%3))), data.frame(c = rep(3, (n%/%3))))

  # Algorithme K-médoïdes maison 
  res_k_medoids_maison <- algo_k_medoids(cluster_test, k)$clusters
  rands$maison[sim] <- (adj.rand.index(cluster_test_real[,], res_k_medoids_maison))

  # Algorithme PAM (K-médoïdes du package cluster)
  res_k_medoids_pam <- pam(cluster_test, k, metric = "euclidean", stand = FALSE, trace.lev = 0)
  rands$pam[sim] <- (adj.rand.index(cluster_test_real[,], res_k_medoids_pam$clustering))

  # Algorithme K-means
  res_k_means <- kmeans(cluster_test, k)
  rands$kmeans[sim] <- (adj.rand.index(cluster_test_real[,], res_k_means$cluster))

  # Algorithme package mclust
  res_mclust <- Mclust(cluster_test, k)
  rands$mclust[sim] <- (adj.rand.index(cluster_test_real[,], res_mclust$classification))

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(rands, aes(x = as.numeric(rownames(rands)))) + 
  geom_point(aes(y = pam, colour = "PAM")) + 
  geom_line(aes(y = pam, colour = "PAM")) + 
  geom_point(aes(y = kmeans, colour = "K-Means")) +
  geom_line(aes(y = kmeans, colour = "K-Means")) +
  geom_point(aes(y = maison, colour = "K-Médoïdes maison")) + 
  geom_line(aes(y = maison, colour = "K-Médoïdes maison")) + 
  geom_point(aes(y = mclust, colour = "Mclust")) +
  geom_line(aes(y = mclust, colour = "Mclust")) +

  xlab("Simulations de test") + ylab("Valeur du Adjusted Rand Index") +
  ggtitle("Evolution du Adjusted Rand Index pour comparer les différentes méthodes")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
df <- iris # Dataset d'origine
df$Species <- as.character(df$Species)
df$Species[df$Species == "setosa"] <- 1
df$Species[df$Species == "versicolor"] <- 2
df$Species[df$Species == "virginica"] <- 3


## ---- results = "hide"--------------------------------------------------------------------------------------------------------------------------------
df_init <- df[1:4] # Dataset sans la dernière colonne (on supprime la dernière colonne)
X <- df_init
k <- 3
cluster_test <- X
cluster_test_real <- as.numeric(df$Species)
library(fossil)
library(cluster)
library(mclust)
n_sim <- 10

rands <- data.frame(matrix(0, nrow = n_sim, ncol = 4))
colnames(rands) <- c("maison", "pam", "kmeans", "mclust")

for (sim in 1:n_sim) {
  # Algorithme K-médoïdes maison TODO
  res_k_medoids_maison <- algo_k_medoids(cluster_test, k)$clusters
  rands$maison[sim] <- (adj.rand.index(cluster_test_real, res_k_medoids_maison))

  # Algorithme PAM (K-médoïdes du package cluster)
  res_k_medoids_pam <- pam(cluster_test, k, metric = "euclidean", stand = FALSE, trace.lev = 0)
  rands$pam[sim] <- (adj.rand.index(cluster_test_real, res_k_medoids_pam$clustering))

  # Algorithme K-means
  res_k_means <- kmeans(cluster_test, k)
  rands$kmeans[sim] <- (adj.rand.index(cluster_test_real, res_k_means$cluster))

  # Algorithme package mclust
  res_mclust <- Mclust(cluster_test, k)
  rands$mclust[sim] <- (adj.rand.index(cluster_test_real, res_mclust$classification))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(rands, aes(x = as.numeric(rownames(rands)))) + 
  geom_point(aes(y = pam, colour = "PAM")) + 
  geom_line(aes(y = pam, colour = "PAM")) + 
  geom_point(aes(y = kmeans, colour = "K-Means")) +
  geom_line(aes(y = kmeans, colour = "K-Means")) +
  geom_point(aes(y = maison, colour = "K-Médoïdes maison")) + 
  geom_line(aes(y = maison, colour = "K-Médoïdes maison")) + 
  geom_point(aes(y = mclust, colour = "Mclust")) +
  geom_line(aes(y = mclust, colour = "Mclust")) +

  xlab("Simulations de test") + ylab("Valeur du Adjusted Rand Index") +
  ggtitle("Evolution du Adjusted Rand Index pour comparer les différentes méthodes")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Importation des bibliothèques externes
library(cluster)
library(factoextra)

# Importation du dataset des iris
df <- iris # Dataset d'origine
df_init <- df[1:4] # Dataset sans la dernière colonne (on supprime la dernière colonne)

# On scale le dataframe
df_init <- scale(df_init)

# Détermination du nombre de clusters avec la méthode des Total Within Sum of Squares
#fviz_nbclust(df_init, pam, method = "wss")

# Le nombre optimal de clusters à appliquer à la méthode PAM est l'abscisse du point où la courbe "se casse"
# Ici, c'est k = 3
k <- 3
# Application de l'algorithme PAM implémentant les k-médoïdes
k_medoids <- pam(df, k, metric = "euclidean", stand = FALSE)

# Visualisation des différentes partitions du dataset sur les premiers plans d'une analyse en composante principale
fviz_cluster(k_medoids, data = df_init, main = "[Algorithme PAM - K-médoïdes] - Cluster plot", outlier.color = "purple")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Détermination du nombre de clusters avec la méthode des Total Within Sum of Squares
#fviz_nbclust(df_init, kmeans, method = "wss")
k <- 3

# Application de l'algorithme kmeans implémentant les kmeans
k_means <- kmeans(df_init, k)

# Visualisation des différentes partitions du dataset sur les premiers plans d'une analyse en composante principale
fviz_cluster(k_means, data = df_init, main = "[Algorithme kmeans - K-means] - Cluster plot", outlier.color = "purple")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
library(mclust)

# Application de l'algorithme mclust (sans spécifier k)
mclust_res <- Mclust(df_init)

# Visualisation des différentes partitions du dataset sur les premiers plans d'une analyse en composante principale
fviz_cluster(mclust_res, data = df_init, main = "[Algorithme Mclust - mclust] - Cluster plot", outlier.color = "purple")

# Application de l'algorithme mclust (avec k = 3)
mclust_res <- Mclust(df_init, 3)

# Visualisation des différentes partitions du dataset sur les premiers plans d'une analyse en composante principale
fviz_cluster(mclust_res, data = df_init, main = "[Algorithme Mclust - mclust] - Cluster plot", outlier.color = "purple")


## -----------------------------------------------------------------------------------------------------------------------------------------------------
# Affichage graphique des clusters
res_k_medoids_maison <- algo_k_medoids(df_init, 3)
fviz_cluster(object = list(data = df_init, cluster = res_k_medoids_maison$clusters), main = "[Algorithme K-médoïdes - implémentation] - Cluster plot", outlier.color = "purple")

