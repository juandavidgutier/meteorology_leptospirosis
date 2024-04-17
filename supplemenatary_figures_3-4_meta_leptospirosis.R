library(dlnm) ; library(mvmeta) ; library(splines) ; library(dplyr); library(ggplot2); library(RCurl)

# LOAD THE DATASET
url_path_leptos = "https://raw.githubusercontent.com/juandavidgutier/meteorology_leptospirosis/master/data_leptos.csv"
muni_Col <- read.csv(url_path_leptos)
dim(muni_Col)
head(muni_Col)

#tertiles
muni_Col$SST12 <- dplyr::ntile(muni_Col$SST12, 3)  
muni_Col$SST3 <- dplyr::ntile(muni_Col$SST3, 3) 
muni_Col$SST34 <- dplyr::ntile(muni_Col$SST34, 3) 
muni_Col$SST4 <- dplyr::ntile(muni_Col$SST4, 3) 
muni_Col$ESOI <- dplyr::ntile(muni_Col$ESOI, 3) 
muni_Col$SOI <- dplyr::ntile(muni_Col$SOI, 3) 
muni_Col$NATL <- dplyr::ntile(muni_Col$NATL, 3) 
muni_Col$SATL <- dplyr::ntile(muni_Col$SATL, 3) 
muni_Col$TROP <- dplyr::ntile(muni_Col$TROP, 3)

# REGIONS
regions <- as.character(unique(muni_Col$municipality)) 

# CREATE A LIST WITH THE REGIONAL SERIES
data <- lapply(regions,function(x) muni_Col[muni_Col$municipality==x,])
names(data) <- regions
m <- length(regions)


# Runoff RANGES
ranges <- t(sapply(data, function(x) range(x$Runoff,na.rm=T)))

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}


################################################################################
#Lag 0-1 months
lag <- c(0,1)

#FIRST STAGE
bound <- colMeans(ranges)

argvar <- list(fun="poly", degree=2, cen=0.01) #df=3
arglag <- list(fun="ns", df=1, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=2, cen=0.01) #df=3
arglag2 <- list(fun="poly", degree=1,intercept=FALSE)
argvar3 <- list(fun="bs", df=2, cen=0.01) #df=3
arglag3 <- list(fun="poly", degree=1, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall3 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep="")))

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Runoff,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "mm" #"mm"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF Runoff AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:100/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF Runoff AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.01,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])



# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Supplementary figure 3A
Runoff <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

Runoff_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(Runoff, best_model, best_model_h, best_model_l))
names <- c("Runoff", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

sf3A = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=Runoff,
                  ymin= best_model_l,  
                  ymax= best_model_h), 
              fill='grey',alpha=0.2) + 
  geom_line(aes(x=Runoff,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black") +
  geom_vline(aes(xintercept=Runoff_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=Runoff_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  ggtitle("A") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf3A)


# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN
(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

#90TH AND 95TH PERCENTILES
print(Runoff_tiles)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["276",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["291",]),3)


#META-REGRESSION 
# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MPI)
MPI <- unique(muni_Col$MPI)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~MPI))
summary(mvallmpi)

#RR for Runoff and MULTIDIMENSIONAL POVERTY INDEX for each city
#exp(mvallmpi[["model"]])

# PREDICTION FROM META-REGRESSION
val <- round(quantile(MPI,c(25,75)/100),1)
predall <- predict(mvallmpi,data.frame(MPI=val),vcov=T)

cpallmpiat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallmpiat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"MPI"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","MPI")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall, Sall, method = "mm") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall ~ MPI, Sall, method = "mm", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

# Supplementary figure 4A
Runoff <- as.data.frame(cpallmpiat75$predvar)
at75 <- as.data.frame(cpallmpiat75$allRRfit)
at75_h <- as.data.frame(cpallmpiat75$allRRhigh)
at75_l <- as.data.frame(cpallmpiat75$allRRlow)
at25 <- as.data.frame(cpallmpiat25$allRRfit)
at25_h <- as.data.frame(cpallmpiat25$allRRhigh)
at25_l <- as.data.frame(cpallmpiat25$allRRlow)
data_runoff <- as.data.frame(cbind(Runoff, at75, at75_h, at75_l, at25, at25_h, at25_l))
names <- c("Runoff", "at75", "at75_h", "at75_l", "at25", "at25_h", "at25_l")
colnames(data_runoff) <-  names

sf4A = ggplot(data_runoff) +
  geom_ribbon(aes(x=Runoff,
                  ymin= at75 - at75_l, #ymin= best_model - conf_int[1],
                  ymax= at75 + at75_h), 
              fill='red',alpha=0.2) + #percentile 90th of MPI
  geom_ribbon(aes(x=Runoff,
                  ymin= at25 - at25_l,
                  ymax= at25 + at25_h),
              fill='blue',alpha=0.2) + #percentile 10th of MPI
  geom_line(aes(x=Runoff,
                y=at75), linewidth=1.05,
            color='red') + #percentile 90th of MPI
  geom_line(aes(x=Runoff,
                y=at25), linewidth=1.05,
            color='blue') + #percentile 10th of MPI
  geom_hline(aes(yintercept=1), color="black") +
  theme_bw() +
  ggtitle("A") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf4A)



####################################################################################
#Lag 0-2 months
lag <- c(0,2)

#FIRST STAGE
bound <- colMeans(ranges)

argvar <- list(fun="poly", degree=2, cen=0.01) #df=3
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.01) #df=3
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.01) #df=3
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
yall3 <- yall2 

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Runoff,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "mm" #"mm"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF Runoff AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:200/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF Runoff AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.01,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])



# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Supplementary figure 3B
Runoff <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

Runoff_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(Runoff, best_model, best_model_h, best_model_l))
names <- c("Runoff", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

sf3B = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=Runoff,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.2) + 
  geom_line(aes(x=Runoff,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black") +
  geom_vline(aes(xintercept=Runoff_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=Runoff_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  ggtitle("B") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf3B)



# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN
(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

#90TH AND 95TH PERCENTILES
print(Runoff_tiles)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["276",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["291",]),3)


#META-REGRESSION 
# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MPI)
MPI <- unique(muni_Col$MPI)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~MPI))
summary(mvallmpi)

#RR for Runoff and MULTIDIMENSIONAL POVERTY INDEX for each city
#exp(mvallmpi[["model"]])

# PREDICTION FROM META-REGRESSION
val <- round(quantile(MPI,c(25,75)/100),1)
predall <- predict(mvallmpi,data.frame(MPI=val),vcov=T)

cpallmpiat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallmpiat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"MPI"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","MPI")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "mm") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall2 ~ MPI, Sall2, method = "mm", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

# Supplementary figure 4B
Runoff <- as.data.frame(cpallmpiat75$predvar)
at75 <- as.data.frame(cpallmpiat75$allRRfit)
at75_h <- as.data.frame(cpallmpiat75$allRRhigh)
at75_l <- as.data.frame(cpallmpiat75$allRRlow)
at25 <- as.data.frame(cpallmpiat25$allRRfit)
at25_h <- as.data.frame(cpallmpiat25$allRRhigh)
at25_l <- as.data.frame(cpallmpiat25$allRRlow)
data_runoff <- as.data.frame(cbind(Runoff, at75, at75_h, at75_l, at25, at25_h, at25_l))
names <- c("Runoff", "at75", "at75_h", "at75_l", "at25", "at25_h", "at25_l")
colnames(data_runoff) <-  names

sf4B = ggplot(data_runoff) +
  geom_ribbon(aes(x=Runoff,
                  ymin= at75 - at75_l, #ymin= best_model - conf_int[1],
                  ymax= at75 + at75_h), 
              fill='red',alpha=0.2) + #percentile 90th of MPI
  geom_ribbon(aes(x=Runoff,
                  ymin= at25 - at25_l,
                  ymax= at25 + at25_h),
              fill='blue',alpha=0.2) + #percentile 10th of MPI
  geom_line(aes(x=Runoff,
                y=at75), linewidth=1.05,
            color='red') + #percentile 90th of MPI
  geom_line(aes(x=Runoff,
                y=at25), linewidth=1.05,
            color='blue') + #percentile 10th of MPI
  geom_hline(aes(yintercept=1), color="black") +
  theme_bw() +
  ggtitle("B") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR")  
print(sf4B)





####################################################################################
#Lag 0-3 months
lag <- c(0,3)

#FIRST STAGE
bound <- colMeans(ranges)

argvar <- list(fun="poly", degree=2, cen=0.01) #df=3
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.01) #df=3
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.01) #df=3
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 

yall3 <- yall2 

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Runoff,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "mm" #"mm"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF Runoff AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:300/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF Runoff AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.01,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])



# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Supplementary figure 3C
Runoff <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

Runoff_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(Runoff, best_model, best_model_h, best_model_l))
names <- c("Runoff", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

sf3C = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=Runoff,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.2) + 
  geom_line(aes(x=Runoff,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black") +
  geom_vline(aes(xintercept=Runoff_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=Runoff_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  ggtitle("C") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR")  
print(sf3C)



# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN
(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

#90TH AND 95TH PERCENTILES
print(Runoff_tiles)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["276",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["291",]),3)


#META-REGRESSION 
# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MPI)
MPI <- unique(muni_Col$MPI)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~MPI))
summary(mvallmpi)

#RR for Runoff and MULTIDIMENSIONAL POVERTY INDEX for each city
#exp(mvallmpi[["model"]])

# PREDICTION FROM META-REGRESSION
val <- round(quantile(MPI,c(25,75)/100),1)
predall <- predict(mvallmpi,data.frame(MPI=val),vcov=T)

cpallmpiat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallmpiat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"MPI"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","MPI")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "mm") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall2 ~ MPI, Sall2, method = "mm", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

# Supplementary figure 4C
Runoff <- as.data.frame(cpallmpiat75$predvar)
at75 <- as.data.frame(cpallmpiat75$allRRfit)
at75_h <- as.data.frame(cpallmpiat75$allRRhigh)
at75_l <- as.data.frame(cpallmpiat75$allRRlow)
at25 <- as.data.frame(cpallmpiat25$allRRfit)
at25_h <- as.data.frame(cpallmpiat25$allRRhigh)
at25_l <- as.data.frame(cpallmpiat25$allRRlow)
data_runoff <- as.data.frame(cbind(Runoff, at75, at75_h, at75_l, at25, at25_h, at25_l))
names <- c("Runoff", "at75", "at75_h", "at75_l", "at25", "at25_h", "at25_l")
colnames(data_runoff) <-  names

sf4C = ggplot(data_runoff) +
  geom_ribbon(aes(x=Runoff,
                  ymin= at75 - at75_l, #ymin= best_model - conf_int[1],
                  ymax= at75 + at75_h), 
              fill='red',alpha=0.2) + #percentile 90th of MPI
  geom_ribbon(aes(x=Runoff,
                  ymin= at25 - at25_l,
                  ymax= at25 + at25_h),
              fill='blue',alpha=0.2) + #percentile 10th of MPI
  geom_line(aes(x=Runoff,
                y=at75), linewidth=1.05,
            color='red') + #percentile 90th of MPI
  geom_line(aes(x=Runoff,
                y=at25), linewidth=1.05,
            color='blue') + #percentile 10th of MPI
  geom_hline(aes(yintercept=1), color="black") +
  theme_bw() +
  ggtitle("C") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf4C)







################################################################################
#Lag 0-4 months
lag <- c(0,4)

#FIRST STAGE
# MAIN MODEL
argvar <- list(fun="poly", degree=2, cen=0.01) #df=3
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.01) #df=3
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.01) #df=3
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 

yall3 <- yall2 


# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Runoff,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})


# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "mm" #"mm"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)



# BASES OF Runoff AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 1:400/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF Runoff AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.01,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])

#print(exp(cpall$coefficients))
#print(exp(confint(cpall, method="Wald")))

#RR for Runoff for each city
#exp(mvall[["model"]])

# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Supplementary figure 3D
Runoff <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

Runoff_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(Runoff, best_model, best_model_h, best_model_l))
names <- c("Runoff", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

sf3D = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=Runoff,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.2) + 
  geom_line(aes(x=Runoff,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black") +
  geom_vline(aes(xintercept=Runoff_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=Runoff_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  ggtitle("D") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf3D)


# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN

(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["276",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["291",]),3)

#META-REGRESSION 

# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MPI)
MPI <- unique(muni_Col$MPI)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~MPI))
summary(mvallmpi)


# PREDICTION FROM META-REGRESSION
val <- round(quantile(MPI,c(25,75)/100),1)
predall <- predict(mvallmpi,data.frame(MPI=val),vcov=T)

cpallmpiat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallmpiat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"MPI"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","MPI")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "mm") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall2 ~ MPI, Sall2, method = "mm", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)


# Supplementary figure 4D
Runoff <- as.data.frame(cpallmpiat75$predvar)
at75 <- as.data.frame(cpallmpiat75$allRRfit)
at75_h <- as.data.frame(cpallmpiat75$allRRhigh)
at75_l <- as.data.frame(cpallmpiat75$allRRlow)
at25 <- as.data.frame(cpallmpiat25$allRRfit)
at25_h <- as.data.frame(cpallmpiat25$allRRhigh)
at25_l <- as.data.frame(cpallmpiat25$allRRlow)
data_runoff <- as.data.frame(cbind(Runoff, at75, at75_h, at75_l, at25, at25_h, at25_l))
names <- c("Runoff", "at75", "at75_h", "at75_l", "at25", "at25_h", "at25_l")
colnames(data_runoff) <-  names

sf4D = ggplot(data_runoff) +
  geom_ribbon(aes(x=Runoff,
                  ymin= at75 - at75_l, #ymin= best_model - conf_int[1],
                  ymax= at75 + at75_h), 
              fill='red',alpha=0.2) + #percentile 90th of MPI
  geom_ribbon(aes(x=Runoff,
                  ymin= at25 - at25_l,
                  ymax= at25 + at25_h),
              fill='blue',alpha=0.2) + #percentile 10th of MPI
  geom_line(aes(x=Runoff,
                y=at75), linewidth=1.05,
            color='red') + #percentile 90th of MPI
  geom_line(aes(x=Runoff,
                y=at25), linewidth=1.05,
            color='blue') + #percentile 10th of MPI
  geom_hline(aes(yintercept=1), color="black") +
  theme_bw() +
  ggtitle("D") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf4D)




################################################################################
#Lag 0-5 months
lag <- c(0,5)

#FIRST STAGE
# MAIN MODEL
argvar <- list(fun="poly", degree=2, cen=0.01) #df=3
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.01) #df=3
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.01) #df=3
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 

yall3 <- yall2 


# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Runoff,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })
    
    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})


# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "mm" #"mm"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)



# BASES OF Runoff AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 1:500/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))

# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF Runoff AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.01,
                      model.link="log",by=0.1,from=bound[1],to=bound[2])

#print(exp(cpall$coefficients))
#print(exp(confint(cpall, method="Wald")))

#RR for Runoff for each city
#exp(mvall[["model"]])

# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Supplementary figure 3E
Runoff <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

Runoff_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(Runoff, best_model, best_model_h, best_model_l))
names <- c("Runoff", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

sf3E = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=Runoff,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.2) + 
  geom_line(aes(x=Runoff,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black") +
  geom_vline(aes(xintercept=Runoff_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=Runoff_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  ggtitle("E") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf3E)


# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN

(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["276",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["291",]),3)

#META-REGRESSION 

# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MPI)
MPI <- unique(muni_Col$MPI)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~MPI))
summary(mvallmpi)


# PREDICTION FROM META-REGRESSION
val <- round(quantile(MPI,c(25,75)/100),1)
predall <- predict(mvallmpi,data.frame(MPI=val),vcov=T)

cpallmpiat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                          model.link="log",by=0.2)
cpallmpiat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                          model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"MPI"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","MPI")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "mm") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall2 ~ MPI, Sall2, method = "mm", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)


# Supplementary figure 4E
Runoff <- as.data.frame(cpallmpiat75$predvar)
at75 <- as.data.frame(cpallmpiat75$allRRfit)
at75_h <- as.data.frame(cpallmpiat75$allRRhigh)
at75_l <- as.data.frame(cpallmpiat75$allRRlow)
at25 <- as.data.frame(cpallmpiat25$allRRfit)
at25_h <- as.data.frame(cpallmpiat25$allRRhigh)
at25_l <- as.data.frame(cpallmpiat25$allRRlow)
data_runoff <- as.data.frame(cbind(Runoff, at75, at75_h, at75_l, at25, at25_h, at25_l))
names <- c("Runoff", "at75", "at75_h", "at75_l", "at25", "at25_h", "at25_l")
colnames(data_runoff) <-  names

sf4E = ggplot(data_runoff) +
  geom_ribbon(aes(x=Runoff,
                  ymin= at75 - at75_l, #ymin= best_model - conf_int[1],
                  ymax= at75 + at75_h), 
              fill='red',alpha=0.2) + #percentile 90th of MPI
  geom_ribbon(aes(x=Runoff,
                  ymin= at25 - at25_l,
                  ymax= at25 + at25_h),
              fill='blue',alpha=0.2) + #percentile 10th of MPI
  geom_line(aes(x=Runoff,
                y=at75), linewidth=1.05,
            color='red') + #percentile 90th of MPI
  geom_line(aes(x=Runoff,
                y=at25), linewidth=1.05,
            color='blue') + #percentile 10th of MPI
  geom_hline(aes(yintercept=1), color="black") +
  theme_bw() +
  ggtitle("E") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf4E)




################################################################################
#Lag 0-6 months
lag <- c(0,6)


argvar <- list(fun="poly", degree=2, cen=0.01) #df=3
arglag <- list(fun="ns", df=2, intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(fun="ns", df=3, cen=0.01) #df=3
arglag2 <- list(fun="poly", degree=2,intercept=FALSE)
argvar3 <- list(fun="bs", df=3, cen=0.01) #df=3
arglag3 <- list(fun="poly", degree=2, intercept=FALSE)

# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES

yall <- matrix(NA,length(data),2,dimnames=list(regions,paste("b",seq(2),sep=""))) 
yall2 <- matrix(NA,length(data),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 

yall3 <- yall2 


# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND boundARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Runoff,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Runoff,lag=lag,argvar=argvar3,arglag=arglag3)
    })
  
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                     family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst2 <- glm(Cases ~ cb2 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    mfirst3 <- glm(Cases ~ cb3 + + Year + Month + SST12 + SST3 + SST34 + SST4 + ESOI + SOI + NATL + SATL + TROP, maxit = 1000,
                      family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
    
    # REDUCTION TO SUMMARY ASSOCIATIONS
    
    # TO OVERALL CUMULATIVE SUMMARY
    suppressWarnings({
      crall <- crossreduce(cb,mfirst)
      crall2 <- crossreduce(cb2,mfirst2)
      crall3 <- crossreduce(cb3,mfirst3)
    })

    # STORE THE RESULTS
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    yall[i,] <- coef(crall)
    Sall[[i]] <- vcov(crall)
    
    # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
    yall2[i,] <- coef(crall2)
    yall3[i,] <- coef(crall3)
    Sall2[[i]] <- vcov(crall2)
    Sall3[[i]] <- vcov(crall3)
    
    # Q-AIC
    qaic[i] <- fqaic(mfirst)
    qaic2[i] <- fqaic(mfirst2)
    qaic3[i] <- fqaic(mfirst3)
    
  }
})


# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "mm" #"mm"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall)

# BASES OF Runoff AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb2,"argvar")))
xlag <- 0:600/100
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2,"arglag")))


# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall2,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF Runoff AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall), cen=0.01,
                   model.link="log",by=0.1,from=bound[1],to=bound[2])


# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Supplementary figure 3F
Runoff <- as.data.frame(cpall$predvar)
best_model <- as.data.frame(cpall$allRRfit)
best_model_h <- as.data.frame(cpall$allRRhigh)
best_model_l <- as.data.frame(cpall$allRRlow)

Runoff_tiles <- round(quantile(cpall$predvar,c(90,95)/100),1)

data_runoff_model <- as.data.frame(cbind(Runoff, best_model, best_model_h, best_model_l))
names <- c("Runoff", "best_model", "best_model_h", "best_model_l")
colnames(data_runoff_model) <-  names

sf3F = ggplot(data_runoff_model) +
  geom_ribbon(aes(x=Runoff,
                  ymin= best_model_l, 
                  ymax= best_model_h), 
              fill='grey',alpha=0.2) + 
  geom_line(aes(x=Runoff,
                y=best_model), linewidth=1.05,
            color='red') + 
  
  geom_hline(aes(yintercept=1), color="black") +
  geom_vline(aes(xintercept=Runoff_tiles[1]), color="dark orange", linetype="dashed") +
  geom_vline(aes(xintercept=Runoff_tiles[2]), color="red", linetype="dashed") +
  theme_bw() +
  ggtitle("F") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf3F)


# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN
(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)

# OVERALL EFFECTS AT TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["276",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["291",]),3)


#META-REGRESSION 

# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MPI)
MPI <- unique(muni_Col$MPI)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~MPI))
summary(mvallmpi)


#RR for Runoff and MULTIDIMENSIONAL POVERTY INDEX for each city
#exp(mvallmpi[["model"]])

# PREDICTION FROM META-REGRESSION
val <- round(quantile(MPI,c(25,75)/100),1)
predall <- predict(mvallmpi,data.frame(MPI=val),vcov=T)

cpallmpiat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
                        model.link="log",by=0.2)
cpallmpiat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
                        model.link="log",by=0.2)

# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"MPI"),3)

## CREATE TABLE FOR ONLY WALD-TEST STATISTICS
tab <- matrix(NA,2,3)
colnames(tab) <- c("stat","df","p")
rownames(tab) <- c("Intercept-only","MPI")

ftab <- function(model,mref=NULL) { ## where, m <- length(datalist)
  if(!is.null(mref)) {
    coef <- coef(model)[-grep("Int",names(coef(model)))]
    vcov <- vcov(model)[-grep("Int",names(coef(model))),-grep("Int",names(coef(model)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
}

mv <- mvmeta(yall2, Sall2, method = "mm") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall2 ~ MPI, Sall2, method = "mm", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)

# Supplementary figure 4F
Runoff <- as.data.frame(cpallmpiat75$predvar)
at75 <- as.data.frame(cpallmpiat75$allRRfit)
at75_h <- as.data.frame(cpallmpiat75$allRRhigh)
at75_l <- as.data.frame(cpallmpiat75$allRRlow)
at25 <- as.data.frame(cpallmpiat25$allRRfit)
at25_h <- as.data.frame(cpallmpiat25$allRRhigh)
at25_l <- as.data.frame(cpallmpiat25$allRRlow)
data_runoff <- as.data.frame(cbind(Runoff, at75, at75_h, at75_l, at25, at25_h, at25_l))
names <- c("Runoff", "at75", "at75_h", "at75_l", "at25", "at25_h", "at25_l")
colnames(data_runoff) <-  names

sf4F = ggplot(data_runoff) +
  geom_ribbon(aes(x=Runoff,
                  ymin= at75 - at75_l, #ymin= best_model - conf_int[1],
                  ymax= at75 + at75_h), 
              fill='red',alpha=0.2) + #percentile 90th of MPI
  geom_ribbon(aes(x=Runoff,
                  ymin= at25 - at25_l,
                  ymax= at25 + at25_h),
              fill='blue',alpha=0.2) + #percentile 10th of MPI
  geom_line(aes(x=Runoff,
                y=at75), linewidth=1.05,
            color='red') + #percentile 90th of MPI
  geom_line(aes(x=Runoff,
                y=at25), linewidth=1.05,
            color='blue') + #percentile 10th of MPI
  geom_hline(aes(yintercept=1), color="black") +
  theme_bw() +
  ggtitle("F") +   # Set the title
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +  # Adjust title size and alignment
  labs(x="Runoff", y="RR") 
print(sf4F)























