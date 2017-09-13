#' @title Creates a Neighbourhood Using Locality Sensitive Hashing for Gaussian Projections
#'
#' @description This package uses Locality Sensitive Hashing and creates a Neighbourhood Graph for a datset and calculates the ARI value for the same. It uses Gaussian Random planes to decide the nature of a given point.
#'
#' @param
#' mydata A data frame consisting of the data set without the class column
#'
#' @param
#' result9 A column which consists of the class column
#'
#' @return NULL
#'
#' @examples LSH_Gaussian(iris[,-5],iris$Species)
#'
#' @export
#' LSH_Gaussian

LSH_Gaussian <- function(mydata,result9)
{
  #Store the dataset without the class column
  mydata1 = mydata
  #Store only the class column of the dataset
  result8 = result9
  result1 <- t(mydata1)
  mode(result1) <- "numeric"
  cols <- c(1:nrow(mydata1))
  rows <- ncol(mydata1)
  set.seed(100)
  x43 <- NULL
  l9 <-list()
  list3 <- list()
  notables <- 10
  list4 <- list()
  l4 <- list()
  hc1<- NULL
  i <- 1
  l5 <- list()
  #------------------------------------------------------------------#

  #LSH Fit Module

  for(x5 in 1:notables)
  {
    d <- floor(runif(1,min = 10, max = 31))
    norm.random.projection <- function(d, mydata1, scaling=TRUE)
    {
      d.original <- nrow(result1)
      # Projection matrix
      P <- rnorm(d*d.original,mean =0,sd = 1);
      P <- matrix(P, nrow=d,ncol = d.original);

      # random data projection
      if (scaling == TRUE)
        reduced.m <- sqrt(1/d) * (P%*%result1)
      else
        reduced.m <- P%*%result1;
      reduced.m
    }
    m2 <- norm.random.projection(d,result1)
    result2 <- sign(m2)
    result2 <- (result2+1)/2
    M2 <- as.matrix(result2)
    M3 <- t(M2)
    hc1<- NULL
    i <- 1
    l5 <- list()
    l4 <- (lapply(seq_len(nrow(M3)), function(i) M3[i,]))
    for(i in 1:nrow(mydata1))
    {
      hc1 <- paste0(l4[[i]],collapse = "")
      l5[[hc1]] = c(l5[[hc1]],i)                    #generating buckets
    }
    list4 <- append(list4,(list(l5)))
  }

  #-----------------------------------------------------------#

  #LSH Query Module
  for(n in 1:nrow(mydata1))
  {
    for(x5 in 1:notables)
    {
      hc2 <- paste0(l4[[n]],collapse = "")
      len3 <- length(names(list4[[x5]]))
      for(j in 1:len3)
      {
        if(hc2 == names(list4[[x5]][j]))
        {
          break()
        }
      }                                             #matching bucket Id
      len4 <- NULL
      len4 <- length(list4[[x5]][[j]])
      disteu <- NULL
      disteu1 <- NULL
      distham <- NULL
      x41 <- NULL
      x42 <- NULL
      k <- 1
      k1 <- 1
      k2 <- 1
      k3 <- 1                                  #declaring variables
      if(len4 > notables)
      {
        for(k in 1:len4)
        {
          #disteu[k] <- dist2(mydata1[n,],mydata1[list4[[x5]][[j]][k],],method = "euclidean")
          disteu[k] <- lsa::cosine(result1[,n],result1[,list4[[x5]][[j]][k]])
        }
        x41 <- order(disteu,na.last = TRUE,decreasing = FALSE)
        x41 <- x41[2:(notables+1)]
        for(k1 in 1:notables)
        {
          x42[k1] <- list4[[x5]][[j]][x41[k1]]
        }
        disteu1 <- sort(disteu,decreasing = FALSE)
        disteu1 <- disteu1[2:(notables+1)]
      } else {
        for(k2 in 1:len3)
        {
          distham[k2] <- stringdist::stringdist(names(list4[[x5]][j]), names(list4[[x5]][k2]), method = c("hamming"))
        }
        disthamor <- order(distham,decreasing = FALSE)
        disthamor <- disthamor[2:len3]
        for(k3 in 1:(len3-1))
        {
          d1 <- disthamor[k3]
          len5 <- length(list4[[x5]][[d1]])
          if(len5 > (notables-1))
          {
            for(k in 1:len5)
            {
              #disteu[k] <- dist2(mydata1[n,],mydata1[list4[[x5]][[d1]][k],],method = "euclidean")

              disteu[k] <- lsa::cosine(result1[,n],result1[,list4[[x5]][[d1]][k]])

            }
            x41 <- order(disteu,na.last = TRUE,decreasing = FALSE)

            for(k1 in 1:notables)
            {
              x42[k1] <- list4[[x5]][[d1]][x41[k1]]
            }
            disteu1 <- sort(disteu,decreasing = FALSE)
            disteu1 <- disteu1[1:notables]
            break()
          }
        }

      }
      x43[[x5]] <- x42
      x54 <- unique(unlist(x43))                     #calculating unique neighbours
      len6 <- length(x54)
      if(len6 >= notables)
        break()
    }
    disteu23 <- NULL
    for(k in 1:len6)
    {
      #disteu23[k] <- dist2(mydata1[n,],mydata1[x54[k],],method = "euclidean")
      disteu23[k] <- lsa::cosine(result1[,n],result1[,x54[k]])
    }
    x45 <- NULL
    x45 <- order(disteu23,na.last = TRUE,decreasing = FALSE)
    x55 <- NULL
    for(k in 1:notables)
    {
      x55[k] <- x54[x45[k]]
    }
    disteu24 <- NULL
    disteu24 <- disteu23[1:notables]
    l9 <- append(l9,list(c(x55)))
  }
  m4 <- NULL
  m4 <- do.call(rbind, l9)                        #generating the neighbourhood matrix
  print("Neighbourhood Matrix")
  print(m4)
  #--------------------------------------------------------#

  #Creating the clusters and drawing neighbourhood graph
  m7 <- melt(m4)
  m7 <- m7[,-2]
  m8 <- as.matrix(m7)
  g1 <- graph_from_edgelist(m8,directed = FALSE)
  plot(g1,vertex.size = 5)
  clstur <- cluster_louvain(g1)
  pred <- rep(0,nrow(mydata1))
  for(y2 in 1:length(clstur))
  {
    pred[clstur[[y2]]] = y2
  }
  randIndex(pred, result8,correct = T,original = F)
}
