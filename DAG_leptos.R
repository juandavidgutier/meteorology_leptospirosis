library(ggdag)
library(dagitty)
library(lavaan)
library(dplyr)
library(GGally)
library(tidyr)



#implied Conditional Independencies
url_path_leptos = "https://raw.githubusercontent.com/juandavidgutier/meteorology_leptospirosis/master/data_leptos.csv"
dataset <- read.csv(url_path_leptos)

dataset <- select(dataset, Cases, Runoff, SST12, SST3, SST34, SST4, ESOI, SOI, NATL, SATL, TROP, Year, Month)
dataset <- dataset[complete.cases(dataset), ] 
str(dataset)


#DAG 
dag <- dagitty('dag {
Cases [pos="0, 0.5"]
Runoff  [pos="-1, 0.5"]

Year [pos="-1.7, -1.2"]
Month [pos="-1.0, -1.0"]

SOI [pos="-1.6, 1.1"]
ESOI [pos="-1.7, 1.2"]
SST3 [pos="-1.8, 1.3"]
SST4 [pos="-1.9, 1.4"]
SST34 [pos="-2, 1.5"]
SST12 [pos="-2.1, 1.6"]
NATL [pos="-2.2, 1.7"]
SATL [pos="-2.3, 1.8"]
TROP [pos="-2.4, 1.9"]


SST12 -> SST3
SST12 -> SST34
SST12 -> SST4
SST12 -> SOI
SST12 -> ESOI
SST12 -> NATL
SST12 -> SATL
SST12 -> TROP

SST3 -> SST34
SST3 -> SST4
SST3 -> SOI
SST3 -> ESOI
SST3 -> NATL
SST3 -> SATL
SST3 -> TROP

SST34 -> SST4
SST34 -> SOI
SST34 -> ESOI
SST34 -> NATL
SST34 -> SATL
SST34 -> TROP

SST4 -> SOI
SST4 -> ESOI
SST4 -> NATL
SST4 -> SATL
SST4 -> TROP

SOI -> ESOI
SOI -> NATL
SOI -> SATL
SOI -> TROP

ESOI -> NATL
ESOI -> SATL
ESOI -> TROP

NATL -> SATL
NATL -> TROP

SATL -> TROP

SST12 -> Runoff
SST3 -> Runoff
SST34 -> Runoff
SST4 -> Runoff
SOI -> Runoff
ESOI -> Runoff
NATL -> Runoff
SATL -> Runoff
TROP -> Runoff

SST12 -> Cases
SST3 -> Cases
SST34 -> Cases
SST4 -> Cases
SOI -> Cases
ESOI -> Cases
NATL -> Cases
SATL -> Cases
TROP -> Cases


Year -> SST12
Year -> SST3
Year -> SST34
Year -> SST4
Year -> SOI
Year -> ESOI
Year -> NATL
Year -> SATL
Year -> TROP
Year -> Runoff
Year -> Cases

Month -> Cases
Month -> Runoff
Month -> SST12
Month -> SST3
Month -> SST34
Month -> SST4
Month -> SOI
Month -> ESOI
Month -> NATL
Month -> SATL
Month -> TROP


Runoff -> Cases

}')  


plot(dag)


## check whether any correlations are perfect (i.e., collinearity)
myCov <- cov(dataset)
round(myCov, 2)

myCor <- cov2cor(myCov)
noDiag <- myCor
diag(noDiag) <- 0
any(noDiag == 1)

## if not, check for multicollinearity (i.e., is one variable a linear combination of 2+ variables?)
det(myCov) < 0
## or
any(eigen(myCov)$values < 0)


## Independencias condicionales
impliedConditionalIndependencies(dag, max.results=3)
corr <- lavCor(dataset)

summary(corr)

#plot
localTests(dag, sample.cov=corr, sample.nobs=nrow(dataset), max.conditioning.variables=3)
plotLocalTestResults(localTests(dag, sample.cov=corr, sample.nobs=nrow(dataset)), xlim=c(-1,1))


#identification
simple_dag <- dagify(
  Cases ~  Runoff + SST12 + SST3 + SST34 + SST4 + SOI + EqSOI + NATL + SATL + TROP  + Year + Month,
  Runoff ~ SST12 + SST3 + SST34 + SST4 + SOI + EqSOI + NATL + SATL + TROP + Year + Month,
  SST12 ~ SST3 + SST34 + SST4 + SOI + EqSOI + NATL + SATL +  TROP + Year + Month,
  SST3 ~ SST34 + SST4 + SOI + EqSOI + NATL + SATL +  TROP + Year + Month,
  SST34 ~ SST4 + SOI + EqSOI + NATL + SATL +  TROP + Year + Month,
  SST4 ~ SOI + EqSOI + NATL + SATL +  TROP + Year + Month,
  SOI ~ EqSOI + NATL + SATL +  TROP + Year + Month,
  EqSOI ~ NATL + SATL +  TROP + Year + Month,
  NATL ~ SATL +  TROP + Year + Month,
  SATL ~  TROP + Year + Month,
  exposure = "Runoff",
  outcome = "Cases",
  coords = list(x = c(Runoff=2, Month=1, Cases=2, SST12=3, SST3=3.1, SST34=3.2, SST4=3.3, SOI=3.4, EqSOI=3.5, NATL=3.6, SATL=3.7, TROP=3.8,
                      Year=3.5),
                y = c(Runoff=2, Month=3, Cases=1, SST12=3, SST3=3.1, SST34=3.2, SST4=3.3, SOI=3.4, EqSOI=3.5, NATL=3.6, SATL=3.7, TROP=3.8,
                      Year=1.8))
)


# theme_dag
ggdag(simple_dag) + 
  theme_dag()

ggdag_status(simple_dag) +
  theme_dag()


#adjust
adjustmentSets(simple_dag,  type = "minimal")

ggdag_adjustment_set(simple_dag, shadow = TRUE) +
  theme_dag()

