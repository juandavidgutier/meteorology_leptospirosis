library(dlnm) ; library(mvmeta) ; library(splines) ; library(dplyr); library(RCurl)

# LOAD THE DATASET
url_path_leptos = "https://raw.githubusercontent.com/juandavidgutier/meteorology_leptospirosis/master/20muni_leptos_09_18.csv"
muni_Col <- read.csv(url_path_leptos)
dim(muni_Col)
head(muni_Col)

# Fig 1
ts_20muni <- muni_Col %>%
        group_by(period) %>%
        summarise(total_cases = sum(Cases))

Cases <- ts(ts_20muni$total_cases, start = c(2009,1), frequency = 52)
plot(Cases)


# REGIONS
regions <- as.character(unique(muni_Col$municipality)) 

# CREATE A LIST WITH THE REGIONAL SERIES
data <- lapply(regions,function(x) muni_Col[muni_Col$municipality==x,])
names(data) <- regions
m <- length(regions)

# RAINFALL RANGES
ranges <- t(sapply(data, function(x) range(x$Rainfall,na.rm=T)))

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}


bound <- colMeans(ranges)
argvar <- list(type="ns",cen=5)
arglag <- list(fun="ns",intercept=FALSE)


each_city <- matrix(NA,20,4)
colnames(each_city) <- c("city","RR","2.5%", "97.5%")

ctr_city = 1

for (i in 1:20){ 
  city_name <- regions[[i]]
  
  each_city[ctr_city,1] = city_name
  
  suppressWarnings(
    cb <- crossbasis(data[[i]]$Rainfall,lag=6,argvar=argvar,arglag=arglag)
  )
  
 
  model <- glm(Cases ~ cb + year + week, 
               family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
  
  
  
  rr <- exp(model$coefficients[2])
  each_city[ctr_city,2] = rr
  
  conf_int <- exp(confint(model, method="Wald"))
  conf_int_i <- conf_int[2,]
  each_city[ctr_city,3] = conf_int_i[1]
  each_city[ctr_city,4] = conf_int_i[2]
  
  ctr_city = ctr_city + 1
  

}

each_city_rain <- as.data.frame(each_city)
each_city_rain$RR <- as.numeric(each_city_rain$RR)
each_city_rain$'2.5%' <- as.numeric(each_city_rain$'2.5%')
each_city_rain$'97.5%' <- as.numeric(each_city_rain$'97.5%')

print(each_city_rain)




#FIRST STAGE
# MAIN MODEL
lag <- c(0,6)
bound <- colMeans(ranges)
argvar <- list(type="ns",cen=5)
arglag <- list(fun="ns",intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(type="ns",cen=5)
arglag2 <- list(fun="bs",intercept=FALSE)
argvar3 <- list(fun="ns",cen=5)
arglag3 <- list(fun="poly", degree=2,intercept=FALSE)


# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES
# 1 is the coefficients in crall, crall2 and crall3, in the original paper of Gasparinni was 4
yall <- matrix(NA,length(data),1,dimnames=list(regions,paste("b",seq(1),sep=""))) 
yall2 <- yall3 <- yall

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Shot <- Scold <- Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0


# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Rainfall,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Rainfall,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Rainfall,lag=lag,argvar=argvar3,arglag=arglag3)
    })
  
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb+year+week,family=quasipoisson(),sub,offset(log(sub$total_population))) # ojo offset
    mfirst2 <- glm(Cases ~ cb2+year+week,family=quasipoisson(),sub,offset(log(sub$total_population)))
    mfirst3 <- glm(Cases ~ cb3+year+week,family=quasipoisson(),sub,offset(log(sub$total_population)))
    
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


# TEST: REDUCTION OF ALTERNATIVE MODELS TO THE SPACE OF THE PREDICTOR RETURNS
suppressWarnings(coef(crosspred(cb3,mfirst3)))
suppressWarnings(coef(crossreduce(cb3,mfirst3)))

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "reml"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall~1,Sall,method=method)
summary(mvall)


print(exp(mvall$coefficients))
print(exp(confint(mvall, method="Wald")))


#RR for rainfall for each city
exp(mvall[["model"]])


# OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
mvall2 <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall2)
mvall3 <- mvmeta(yall3~1,Sall3,method=method)
summary(mvall3)


# CREATE BASES FOR PREDICTION

# BASES OF RAINFALL AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb,"argvar")))
xlag <- 0:300/50
blag <- do.call("onebasis",c(list(x=xlag),attr(cb,"arglag")))


# REGION-SPECIFIC FIRST-STAGE SUMMARIES
regall <- apply(yall,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF RAINFALL AND LAG VALUES
# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                   model.link="log",by=0.1,from=bound[1],to=bound[2])


print(exp(cpall$coefficients))
print(exp(confint(cpall, method="Wald")))


#RR for rainfall for each city
exp(mvall[["model"]])


# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR ALTERNATIVE MODELS
cpall2 <- crosspred(bvar,coef=coef(mvall2),vcov=vcov(mvall2),
                    model.link="log",by=0.1,from=bound[1],to=bound[2])
cpall3 <- crosspred(bvar,coef=coef(mvall3),vcov=vcov(mvall3),
                    model.link="log",by=0.1,from=bound[1],to=bound[2])


# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Fig 2A
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
plot(cpall,type="n",ylab="RR",xlab="Rainfall (mm)") 
abline(h=1)
lines(cpall,col=2,lwd=2.5)
plot(cpall,ylab="RR",col=2,lwd=2.5,xlab="Rainfall (mm)") 
lines(cpall2,col=3,lty=2,lwd=2.5)
lines(cpall3,col=4,lty=4,lwd=2.5)
legend ("top",c("N-spline (with 95%CI)","B-spline",
                "Poly dg=2"),lty=c(1,2,4),lwd=1.5,col=2:4,bty="n",inset=0.05,
        cex=0.8)



# POINT OF MINIMUM INCIDENCE
cpall$predvar[which.min(cpall$allRRfit)]
round(sum(muni_Col$Rainfall<3.8)/nrow(muni_Col)*100,1)

# Q TEST AND I-SQUARE
# I-SQUARE > 75 SIGNIFICA ALTA HETEROGENEIDAD Y NO DEBERIA HABERSE HECHO EL META-ANALISIS, EN ESTE CASO QUE I-SQUARE=35.3 ESTA BIEN

(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)
(qall2 <- qtest(mvall2))
round(((qall2$Q-qall2$df)/qall2$Q)[1]*100,1)
(qall3 <- qtest(mvall3))
round(((qall3$Q-qall3$df)/qall3$Q)[1]*100,1)


#META-REGRESSION 

# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MP.index)
mp.index <- unique(muni_Col$MP.index)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~mp.index))
summary(mvallmpi)


print(exp(mvallmpi$coefficients))
print(exp(confint(mvallmpi, method="Wald")))


#RR for rainfall and MULTIDIMENSIONAL POVERTY INDEX for each city
exp(mvallmpi[["model"]])


# PREDICTION FROM META-REGRESSION
val <- round(quantile(mp.index,c(10,90)/100),1)
predall <- predict(mvallmpi,data.frame(mp.index=val),vcov=T)

cpallmpiat10 <- crosspred(bvar,coef=predall[1],vcov=matrix(predall[3]),
                        model.link="log",by=0.2)
cpallmpiat90 <- crosspred(bvar,coef=predall[2],vcov=matrix(predall[4]),
                        model.link="log",by=0.2)



# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qallmpi <- qtest(mvallmpi))
#round(((qallmpi$Q-qallmpi$df)/qallmpi$Q)[1]*100,1)


# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvallmpi,"mp.index"),3)


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

mv <- mvmeta(yall, Sall, method = "reml") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall ~ mp.index, Sall, method = "reml", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)





# Fig 3A 
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
#layout(matrix(c(0,1,1,0,2,2,3,3),2,4,byrow=TRUE))

plot(cpallmpiat90,type="n",ci="n",ylab="RR",xlab="Rainfall (mm)")
lines(cpallmpiat10,col=4,lty=1,lwd=2.5,ci="area",
      ci.arg=list(density=20,col=4))
lines(cpallmpiat90,col=2,lty=1,lwd=2.5,ci="area",
      ci.arg=list(density=20,angle=-45,col=2))
legend("top",paste(val),lty=1,col=c(4,2),lwd=2.5,bty="n",
       title="Percentile 10th and 90th of MPI",inset=0.1,cex=0.9)
#mtext("Overall cumulative summary",cex=0.7)






















#################################################################################
#################################################################################
# TEMPERATURE RANGES
ranges <- t(sapply(data, function(x) range(x$Temperature,na.rm=T)))

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}



# ARGUMENTS AND LISTS FOR CROSS-BASIS DEFINITION
bound <- colMeans(ranges)
argvar <- list(type="ns",cen=10)
arglag <- list(fun="ns",intercept=FALSE)

each_city <- matrix(NA,20,4)
colnames(each_city) <- c("city","RR","2.5%", "97.5%")

ctr_city = 1

for (i in 1:20){ 
  city_name <- regions[[i]]
  
  each_city[ctr_city,1] = city_name
  
  suppressWarnings(
    cb <- crossbasis(data[[i]]$Temperature,lag=6,argvar=argvar,arglag=arglag)
  )
  
  
  model <- glm(Cases ~ cb + year + week, 
               family=quasipoisson(), data[[i]], offset = log(data[[i]]$total_population))
  
  
  
  rr <- exp(model$coefficients[2])
  each_city[ctr_city,2] = rr
  
  conf_int <- exp(confint(model, method="Wald"))
  conf_int_i <- conf_int[2,]
  each_city[ctr_city,3] = conf_int_i[1]
  each_city[ctr_city,4] = conf_int_i[2]
  
  ctr_city = ctr_city + 1
  
  
}

each_city_temp <- as.data.frame(each_city)
each_city_temp$RR <- as.numeric(each_city_temp$RR)
each_city_temp$'2.5%' <- as.numeric(each_city_temp$'2.5%')
each_city_temp$'97.5%' <- as.numeric(each_city_temp$'97.5%')

print(each_city_temp)


#FIRST STAGE
# MAIN MODEL
lag <- c(0,6)
bound <- colMeans(ranges)
argvar <- list(type="ns",cen=10)
arglag <- list(fun="ns",intercept=FALSE)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE BUT DIFFERENT LAG SPACE
argvar2 <- list(type="ns",cen=10)
arglag2 <- list(fun="bs",intercept=FALSE)
argvar3 <- list(fun="ns",cen=10)
arglag3 <- list(fun="poly", degree=2,intercept=FALSE)


# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES
# 1 is the coefficients in crall, crall2 and crall3, in the original paper of Gasparinni was 4
yall <- matrix(NA,length(data),1,dimnames=list(regions,paste("b",seq(1),sep=""))) 
yall2 <- yall3 <- yall

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Shot <- Scold <- Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0


# RUN THE MODEL FOR EACH CITY

# LOOP FOR CITIES
# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
system.time({
  for(i in seq(data)) {
    
    # LOAD
    sub <- data[[i]]
    
    # DEFINE THE CROSS-BASES
    suppressWarnings({
      cb <- crossbasis(sub$Temperature,lag=lag,argvar=argvar,arglag=arglag)
      cb2 <- crossbasis(sub$Temperature,lag=lag,argvar=argvar2,arglag=arglag2)
      cb3 <- crossbasis(sub$Temperature,lag=lag,argvar=argvar3,arglag=arglag3)
    })
    
    # RUN THE FIRST-STAGE MODELS
    mfirst <- glm(Cases ~ cb+year+week,family=quasipoisson(),sub,offset(log(sub$total_population)))
    mfirst2 <- glm(Cases ~ cb2+year+week,family=quasipoisson(),sub,offset(log(sub$total_population)))
    mfirst3 <- glm(Cases ~ cb3+year+week,family=quasipoisson(),sub,offset(log(sub$total_population)))
    
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


# TEST: REDUCTION OF ALTERNATIVE MODELS TO THE SPACE OF THE PREDICTOR RETURNS
suppressWarnings(coef(crosspred(cb2,mfirst2)))
suppressWarnings(coef(crossreduce(cb3,mfirst3)))

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)


#SECOND STAGE
# PERFORM MULTIVARIATE META-ANALYSIS

# SELECT THE ESTIMATION METHOD
method <- "reml"

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall~1,Sall,method=method)
summary(mvall)


print(exp(mvall$coefficients))
print(exp(confint(mvall, method="Wald")))


#RR for temperature or each city
exp(mvall[["model"]])


# OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
mvall2 <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall2)
mvall3 <- mvmeta(yall3~1,Sall3,method=method)
summary(mvall3)


# CREATE BASES FOR PREDICTION

# BASES OF TEMPERATURE AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)

bvar <- do.call("onebasis",c(list(x=xvar),attr(cb,"argvar")))
xlag <- 0:300/50
blag <- do.call("onebasis",c(list(x=xlag),attr(cb,"arglag")))


# REGION-SPECIFIC FIRST-STAGE SUMMARIES

regall <- apply(yall,1,function(x) exp(bvar%*%x))

# PREDICTION FOR A GRID OF TEMPERATURE AND LAG VALUES

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                   model.link="log",by=0.1,from=bound[1],to=bound[2])


print(exp(cpall$coefficients))
print(exp(confint(cpall, method="Wald")))


#RR for temperature for each city
exp(mvall[["model"]])


# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR ALTERNATIVE MODELS
cpall2 <- crosspred(bvar,coef=coef(mvall2),vcov=vcov(mvall2),
                    model.link="log",by=0.1,from=bound[1],to=bound[2])
cpall3 <- crosspred(bvar,coef=coef(mvall3),vcov=vcov(mvall3),
                    model.link="log",by=0.1,from=bound[1],to=bound[2])


# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# Fig 2B
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
plot(cpall,type="n",ylab="RR",xlab="Temperature (°C)") 
abline(h=1)
lines(cpall,col=2,lwd=2.5)
plot(cpall,ylab="RR",col=2,lwd=2.5,xlab="Temperature (°C)") 
lines(cpall2,col=3,lty=2,lwd=2.5)
lines(cpall3,col=4,lty=4,lwd=2.5)
legend ("top",c("N-spline (with 95%CI)","B-spline",
                "Poly dg=2"),lty=c(1,2,4),lwd=1.5,col=2:4,bty="n",inset=0.05,
        cex=0.8)



# POINT OF MINIMUM INCIDENCE
cpall$predvar[which.min(cpall$allRRfit)]
round(sum(muni_Col$Temperature<3.8)/nrow(muni_Col)*100,1)

# Q TEST AND I-SQUARE
(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)
(qall2 <- qtest(mvall2))
round(((qall2$Q-qall2$df)/qall2$Q)[1]*100,1)
(qall3 <- qtest(mvall3))
round(((qall3$Q-qall3$df)/qall3$Q)[1]*100,1)

#META-REGRESSION

# INPUT THE META-VARIABLE: MULTIDIMENSIONAL POVERTY INDEX (MP.index)
mp.index <- unique(muni_Col$MP.index)

# MULTIVARIATE META-REGRESSION
(mvallmpi <- update(mvall,.~mp.index))
summary(mvallmpi)


print(exp(mvallmpi$coefficients))
print(exp(confint(mvallmpi, method="Wald")))


#RR for temperature and MULTIDIMENSIONAL POVERTY INDEX for each city
exp(mvallmpi[["model"]])


# PREDICTION FROM META-REGRESSION
val <- round(quantile(mp.index,c(10,90)/100),1)
predall <- predict(mvallmpi,data.frame(mp.index=val),vcov=T)

cpallmpiat10 <- crosspred(bvar,coef=predall[1],vcov=matrix(predall[3]),
                          model.link="log",by=0.2)
cpallmpiat90 <- crosspred(bvar,coef=predall[2],vcov=matrix(predall[4]),
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
round(fwald(mvallmpi,"mp.index"),3)


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

mv <- mvmeta(yall, Sall, method = "reml") ## NO META-PREDICTOR
mvmpi <- mvmeta(yall ~ mp.index, Sall, method = "reml", ) ## INCLUSION OF MULTIDIMENTIONAL POVERTY INDEX AS META-PREDICTOR

## FILL VALUE IN THE TABLE CREATED
tab[1,] <- ftab(mv)
tab[2,] <- ftab(mvmpi, mv) ## THIS WILL AUTOMATICALLY ADD WALD-TEST STATISTIC, DEGREE OF FREEDOM, AND P-VALUE INTO THE TABLE
print(tab)


# Fig 3B 
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
#layout(matrix(c(0,1,1,0,2,2,3,3),2,4,byrow=TRUE))

plot(cpallmpiat90,type="n",ci="n",ylab="RR",xlab="Temperature (°C)")
lines(cpallmpiat10,col=4,lty=1,lwd=2.5,ci="area",
      ci.arg=list(density=20,col=4))
lines(cpallmpiat90,col=2,lty=1,lwd=2.5,ci="area",
      ci.arg=list(density=20,angle=-45,col=2))
legend("top",paste(val),lty=1,col=c(4,2),lwd=2.5,bty="n",
       title="Percentile 10th and 90th of MPI",inset=0.1,cex=0.9)



