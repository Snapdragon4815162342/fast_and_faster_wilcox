library('microbenchmark')
library('bigmemory')
library('gficf')
library('peakRAM')
library('ggplot2')

library('Rcpp')
sourceCpp("Rcpp_FasterWilcox.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp("Rcpp_FastWilcox.cpp", verbose = TRUE, rebuild = TRUE)

# Test to see if the functions work
v1<- c(45, 33, 35, 39, 42)
v2<- c(34, 36, 41, 43, 44, 37)

a <- wilcox.test(v1, v2, alternative = "two.sided", paired = F, exact = F, correct = F)
b <- fasterWilcox(v1, v2, verbose = T)
c <- fastWilcox(v1, v2, verbose = T)

a$p.value  #0.855
b[1]       #0.855
c[1]       #0.855

rm(v1, v2, a, b, c)

set.seed(123)

ndecadi = 4
nRepetitions = 100

# Sul mio computer(usando una virtual machine), questo ciclo for impiega circa  3 minuti

exec_time_gficf <- rep(0, ndecadi)
sd_time_gficf   <- rep(0, ndecadi)
peak_ram_gficf  <- rep(0, ndecadi)

exec_time_fastWilcox <- rep(0, ndecadi)
sd_time_fastWilcox   <- rep(0, ndecadi)
peak_ram_fastWilcox  <- rep(0, ndecadi)

exec_time_fasterWilcox <- rep(0, ndecadi)
sd_time_fasterWilcox   <- rep(0, ndecadi)
peak_ram_fasterWilcox  <- rep(0, ndecadi)



for (i in 1:ndecadi){
  ncol = 10^i;
  # Inizializzo due vettori con valori random
  #M <- big.matrix(nrow_max, ncol_max, type = 'double', init = 0)
  
  a <- rnorm(ncol, mean = 1, sd = 1)
  b <- rnorm(ncol, mean = 1, sd = 1)
  
  # Eseguo fastWilcox, fasterWilcox(le funzioni scritta da me) e rcpp_WMU_test(da gficf) tra i 2 vettori a e b salvando il tempo di esecuzione
  
  #gficf rcpp_WMU_test
  m <- t(as.matrix(c(a,b)))
  
  mb_res_gficf <- microbenchmark(gficf:::rcpp_WMU_test(m,1:length(a),(length(a)+1):length(m)), times = nRepetitions)
  pr_gficf <- peakRAM(gficf:::rcpp_WMU_test(m,1:length(a),(length(a)+1):length(m)))
  
  exec_time_gficf[i] <- mean(mb_res_gficf[,2])/10^9
  sd_time_gficf[i]   <- sd(mb_res_gficf[,2])/10^9
  peak_ram_gficf[i]  <- pr_gficf[1,4]
  
  
  #rcpp_fastWilcox
  
  mb_res_fastWilcox <- microbenchmark(fastWilcox(a, b), times = nRepetitions)
  pr_fastWilcox <- peakRAM(fastWilcox(a, b))
  
  exec_time_fastWilcox[i] <- mean(mb_res_fastWilcox[,2])/10^9
  sd_time_fastWilcox[i] <- sd(mb_res_fastWilcox[,2])/10^9
  peak_ram_fastWilcox[i] <- pr_fastWilcox[1,4] 
  
  
  #rcpp_fasterWilcox
  
  mb_res_fasterWilcox <- microbenchmark(fasterWilcox(a, b), times = nRepetitions)
  pr_fasterWilcox <- peakRAM(fasterWilcox(a, b))
  
  exec_time_fasterWilcox[i] <- mean(mb_res_fasterWilcox[,2])/10^9
  sd_time_fasterWilcox[i] <- sd(mb_res_fasterWilcox[,2])/10^9
  peak_ram_fasterWilcox[i] <- pr_fasterWilcox[1,4] 
  
  print(c("vector size 10^", toString(i), " done..."))  
  
}

# Rimuovo a,b e m perchè occupano troppa memoria
rm(a)
rm(b)
rm(m)


# Eseguo rcpp_parallel_WMU(...) tra le 2 matrici m1 e m2 e salvo il tempo di esecuzione

exec_time_paral   <-rep(0, ndecadi)
sd_time_paral     <-rep(0, ndecadi)
peak_ram_paral    <-rep(0, ndecadi)


for (i in 1:ndecadi){
  m1 <- matrix(1:10, nrow = 10^i, ncol = 10)
  m2 <- matrix(11:20, nrow = 10^i, ncol = 10)
  mb_res_paral <- microbenchmark(gficf:::rcpp_parallel_WMU_test(m1,
                                                                m2,
                                                                printOutput = FALSE),
                                 times = nRepetitions )
  pr_res_paral <- peakRAM(gficf:::rcpp_parallel_WMU_test(m1,
                                                         m2,
                                                         printOutput = FALSE))
  

  exec_time_paral[i] = mean(mb_res_paral[,2])/10^9
  sd_time_paral[i]   = sd(mb_res_paral[,2])/10^9
  peak_ram_paral[i]  = pr_res_paral[1,4]
  
  print(c("vector size 10^", toString(i), " done..."))  
}



# Creo il data.frame
vector_size = 10^c(1:ndecadi)

results <- data.frame(vector_size, 
                      exec_time_gficf, sd_time_gficf, peak_ram_gficf,
                      exec_time_paral, sd_time_paral, peak_ram_paral,
                      exec_time_fasterWilcox, sd_time_fasterWilcox, peak_ram_fasterWilcox,
                      exec_time_fastWilcox, sd_time_fastWilcox, peak_ram_fastWilcox)
names(results) <- c('vector_size', 
                    'time_gficf', 'sd_gficf', 'peak_ram_gficf',
                    'time_paral', 'sd_paral', 'peak_ram_paral',
                    'time_fasterWilcox', 'sd_fasterWilcox', 'peak_ram_fasterWilcox',
                    'time_fastWilcox', 'sd_fastWilcox', 'peak_ram_fastWilcox')

write.csv(results, "benchmark_MWU_test.csv", row.names=FALSE)


# Visualizzo i tempi di esecuzione

ggplot(data = results,
       mapping = aes(x = log10(vector_size)))+
  
  geom_point(aes(y = exec_time_gficf),
             size = 3,
             colour = 'green')+
  geom_line(y = exec_time_gficf, colour= 'green')+
  geom_errorbar( aes(ymin = exec_time_gficf-sd_time_gficf, 
                     ymax = exec_time_gficf+sd_time_gficf), 
                 width=0.4, 
                 colour="green", 
                 alpha=0.9, 
                 size=1.3)+
  
  geom_point(aes(y = exec_time_paral),
             size = 3,
             colour = 'green4')+
  geom_line(y = exec_time_paral, colour= 'green4')+
  geom_errorbar( aes(ymin = exec_time_paral-sd_time_paral, 
                     ymax = exec_time_paral+sd_time_paral), 
                 width=0.4, 
                 colour="green4", 
                 alpha=0.9, 
                 size=1.3)+
  geom_point(aes(y = exec_time_fastWilcox),
             size = 3,
             colour = 'cyan')+
  geom_line(y = exec_time_fastWilcox, colour= 'cyan')+
  geom_errorbar( aes(ymin = exec_time_fastWilcox-sd_time_fastWilcox, 
                     ymax = exec_time_fastWilcox+sd_time_fastWilcox), 
                 width=0.4, 
                 colour="cyan", 
                 alpha=0.9, 
                 size=1.3)+
  
  geom_point(aes(y = exec_time_fasterWilcox),
             size = 3,
             colour = 'blue')+
  geom_line(y = exec_time_fasterWilcox, colour= 'blue')+
  geom_errorbar( aes(ymin = exec_time_fasterWilcox-sd_time_fasterWilcox, 
                     ymax = exec_time_fasterWilcox+sd_time_fasterWilcox), 
                 width=0.4, 
                 colour="blue", 
                 alpha=0.9, 
                 size=1.3)+
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=14, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=14, angle=45))+
  scale_x_continuous(labels = as.character(log10(vector_size)), breaks = log10(vector_size))+
  labs(x = 'log10(input vector size)', y = 'execution time [s]')+
  labs(title = "Execution time (green = gficf, green4 = gficf_parallel, cyan = fasterWilcox)")


# Controllo che le funzioni fastWilcox, fasterWilcox e wilcox restituiscano lo stesso risultato

nDifferentVectorLengths <- 5

numRow = 1000

rmse_fast <- rep(0,nDifferentVectorLengths)
rmse_faster <- rep(0,nDifferentVectorLengths)
rmse_gficf <- rep(0,nDifferentVectorLengths)

pValuesWilcox <- matrix(0, nrow = numRow, ncol = nDifferentVectorLengths)
pValuesGficf <- matrix(0, nrow = numRow, ncol = nDifferentVectorLengths)
pValuesfasterWilcox <- matrix(0, nrow = numRow, ncol = nDifferentVectorLengths)
pValuesfastWilcox   <- matrix(0, nrow = numRow, ncol = nDifferentVectorLengths)



for (k in 1:nDifferentVectorLengths){
  
  numCol = 10^k
  
  randomMat <- matrix(runif(numRow*numCol), nrow = numRow, ncol = numCol)  
  # Anche usando rnorm e rbinom, i risultati sono molto simili (stessa forma dei grafici e stessi ordini di grandezza)
  
  for(i in 1:numRow){
    temp <- wilcox.test(randomMat[i,1:floor(numCol/2)], randomMat[i,(floor(numCol/2)+1):numCol], alternative = "two.sided", paired = F, exact = F, correct = F)
    pValuesWilcox[i,k] <- temp$p.value
    
    
    temp <- gficf:::rcpp_WMU_test(t(as.matrix(randomMat[i,])),1:floor(numCol/2),(floor(numCol/2)+1):numCol)
    pValuesGficf[i,k] <- temp[1]
    
    
    temp <- fastWilcox(randomMat[i,1:floor(numCol/2)], randomMat[i,(floor(numCol/2)+1):numCol])
    pValuesfastWilcox[i,k] <- temp[1]
    
    temp <- fasterWilcox(randomMat[i,1:floor(numCol/2)], randomMat[i,(floor(numCol/2)+1):numCol])
    pValuesfasterWilcox[i,k] <- temp[1]
  }
  
  rmse_fast[k]   <- sqrt(mean((pValuesWilcox[,k] - pValuesfastWilcox[,k])/numRow)^2)
  rmse_faster[k] <- sqrt(mean((pValuesWilcox[,k] - pValuesfasterWilcox[,k])/numRow)^2)
  rmse_gficf[k]  <- sqrt(mean((pValuesWilcox[,k] - pValuesGficf[,k])/numRow)^2)

  print(c("vector size 10^", toString(k), " done... repeated ", toString(numRow), "times"))  
}


library("GMCM")
#, ylim = c(0.4, 0.6)
plot(colMeans(pValuesWilcox))
lines(colMeans(pValuesWilcox))
arrows(x0 = 1:nDifferentVectorLengths,
       y0=colMeans(pValuesWilcox)-GMCM:::colSds(pValuesWilcox),
       x1 = 1:nDifferentVectorLengths,
       y1=colMeans(pValuesWilcox)+GMCM:::colSds(pValuesWilcox),
       code=3,
       angle=90,
       length=0.1)

points(colMeans(pValuesfastWilcox), col= "cyan")
lines(colMeans(pValuesfastWilcox), col= "cyan")
arrows(x0 = 1:nDifferentVectorLengths,
       y0=colMeans(pValuesfastWilcox)-GMCM:::colSds(pValuesfastWilcox),
       x1 = 1:nDifferentVectorLengths,
       y1=colMeans(pValuesfastWilcox)+GMCM:::colSds(pValuesfastWilcox),
       code=3,
       angle=90,
       length=0.1,
       col = "cyan")

points(colMeans(pValuesfasterWilcox), col= "blue")
lines(colMeans(pValuesfasterWilcox), col= "blue")
arrows(x0 = 1:nDifferentVectorLengths,
       y0=colMeans(pValuesfasterWilcox)-GMCM:::colSds(pValuesfasterWilcox),
       x1 = 1:nDifferentVectorLengths,
       y1=colMeans(pValuesfasterWilcox)+GMCM:::colSds(pValuesfasterWilcox),
       code=3,
       angle=90,
       length=0.1,
       col = "blue")

points(colMeans(pValuesGficf), col= "green")
lines(colMeans(pValuesGficf), col= "green")
arrows(x0 = 1:nDifferentVectorLengths,
       y0=colMeans(pValuesGficf)-GMCM:::colSds(pValuesGficf),
       x1 = 1:nDifferentVectorLengths,
       y1=colMeans(pValuesGficf)+GMCM:::colSds(pValuesGficf),
       code=3,
       angle=90,
       length=0.1,
       col = "green")


#ylim = c(-10e-5, 10e-5)

plot(rmse_faster, col = "blue")
lines(rmse_faster, col = "blue")

points(rmse_fast, col = "cyan")
lines(rmse_fast, col = "cyan")

points(rmse_fast, col = "green")
lines(rmse_fast, col = "green")




# Error over p value
plot(1:nrow(pValuesfasterWilcox), (sort(pValuesWilcox[,5]) - sort(pValuesfasterWilcox[,5])),
     #ylim = c(-0.01, 0.01),
     main = paste('Differenza tra Pwilcox e PfasterWilcox al variare del pValue, vectLen = 10^',toString(k)),
     xlab = 'Pvalue[ci sono 10000 campioni, x=0-> pVal=0, x=10000-> pVal=1]',
     ylab = 'Pwilcox - PfasterWilcox')  # vectLen = (10^k)/2

plot(1:nrow(pValuesfastWilcox), (sort(pValuesWilcox[,5]) - sort(pValuesfastWilcox[,5])),
     #ylim = c(-0.01, 0.01),
     main = paste('Differenza tra Pwilcox e PfasterWilcox al variare del pValue, vectLen = 10^',toString(k)),
     xlab = 'Pvalue[ci sono 10000 campioni, x=0-> pVal=0, x=10000-> pVal=1]',
     ylab = 'Pwilcox - PfasterWilcox')  # vectLen = (10^k)/2

plot(1:nrow(pValuesGficf), (sort(pValuesWilcox[,5]) - sort(pValuesGficf[,5])),
     #ylim = c(-0.01, 0.01),
     main = paste('Differenza tra Pwilcox e PfasterWilcox al variare del pValue, vectLen = 10^',toString(k)),
     xlab = 'Pvalue[ci sono 10000 campioni, x=0-> pVal=0, x=10000-> pVal=1]',
     ylab = 'Pwilcox - PfasterWilcox')  # vectLen = (10^k)/2


# Trying to use a big.matrix

bm <- matrix(rnorm(100), nrow = 10, ncol = 10)

bm <- as.big.matrix(bm)

fasterWilcox(bm[3,1:5], bm[3,6:10])

