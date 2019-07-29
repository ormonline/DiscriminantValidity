#+eval=FALSE
library(lavaan)
library(semTools)
library(psych)

# Set up the covariance data

ECSICovData <- matrix(NA,11,11)
rownames(ECSICovData) <- 
  colnames(ECSICovData) <- c("ACSI1", "ACSI2", "ACSI3", "CUEX1", "CUEX2", 
                             "CUEX3", "PERQ1", "PERQ2", "PERQ3", "PERV1",
                             "PERV2")

ECSICovData[upper.tri(ECSICovData, diag = TRUE)] <-
  c(4.00,
    3.23, 4.41,
    2.66, 2.67, 3.61,
    1.81, 1.50, 1.56, 4.41,
    1.85, 1.62, 1.63, 2.63, 4.84,
    1.24, 1.16, 1.09, 1.55, 1.72, 5.29,
    3.04, 2.83, 2.35, 2.07, 1.96, 1.31, 3.61,
    2.81, 2.57, 2.19, 1.55, 1.74, 1.12, 2.67, 3.24,
    1.73, 1.64, 1.32, 0.89, 1.05, 1.41, 1.62, 1.62, 2.89,
    2.66, 2.49, 2.12, 1.40, 1.43, 0.95, 2.22, 2.01, 1.25, 3.24,
    2.99, 2.86, 2.42, 1.52, 1.50, 1.01, 2.34, 2.14, 1.31, 3.05, 4.84)

ECSICovData[lower.tri(ECSICovData)] <-
  t(ECSICovData)[lower.tri(ECSICovData)]

N <- 10417


##### Scale score correlations #######

# A matrix that defines how the indicators are summed

W <- diag(4)[c(1,1,1,2,2,2,3,3,3,4,4),]
rownames(W) <- rownames(ECSICovData)
colnames(W) <- c("ACSI","CUEX","PERQ","PERV")

scaleScoreCorrelations <- cov2cor(t(W)%*%ECSICovData%*%W)

scaleScoreCorrelationsStd <- cov2cor(t(W)%*%cov2cor(ECSICovData)%*%W)

# Scale score correlatins can also be calculated using unit weights,
# which can be implemented by standardizing the indicators before taking the
# sum

##### Disattenuated correlation using tau-equivalent reliability #####

tauEquivRel <- c(alpha(ECSICovData[1:3,1:3])$total$raw_alpha,
                 alpha(ECSICovData[4:6,4:6])$total$raw_alpha,
                 alpha(ECSICovData[7:9,7:9])$total$raw_alpha,
                 alpha(ECSICovData[10:11,10:11])$total$raw_alpha)

disattenuatedCorrelations <- scaleScoreCorrelations/
  sqrt(tauEquivRel%o%tauEquivRel)

diag(disattenuatedCorrelations) <- 1

disattenuatedCorrelations

##### Disattenuated correlation using parallel reliability (HTMT) #####

parallelRel <- c(alpha(ECSICovData[1:3,1:3])$total$std.alpha,
                 alpha(ECSICovData[4:6,4:6])$total$std.alpha,
                 alpha(ECSICovData[7:9,7:9])$total$std.alpha,
                 alpha(ECSICovData[10:11,10:11])$total$std.alpha)

# We use the sum of standardized indicators

disattenuatedCorrelations <- scaleScoreCorrelationsStd/
  sqrt(parallelRel%o%parallelRel)

diag(disattenuatedCorrelations) <- 1

disattenuatedCorrelations

# Alternatively, using semTools
CFAmodel <- "ACSI =~ ACSI1 + ACSI2 + ACSI3
CUEX =~ CUEX1 + CUEX2 + CUEX3
PERQ =~ PERQ1 + PERQ2 + PERQ3
PERV =~ PERV1 + PERV2"

htmt(CFAmodel, sample.cov = ECSICovData)


##### Disattenuated correlation using congeneric reliability #####

CFAest <- cfa(CFAmodel, sample.cov = ECSICovData, sample.nobs = N,
              std.lv = TRUE) # use alternative scaling)

# Use the semTools reliability function
congenericRel <- reliability(CFAest)["omega",1:4]

disattenuatedCorrelations <- scaleScoreCorrelations/
  sqrt(congenericRel%o%congenericRel)

diag(disattenuatedCorrelations) <- 1

disattenuatedCorrelations


##### Cross-loadings (obtained from exploratory factor analysis) #####


# Two factor solution, varimax rotation

EFAest <- fa(ECSICovData, nfactors = 2, rotate = "varimax")

# Pattern coefficients are printed out by default
EFAest

# Structure coefficients equal  pattern coefficients
EFAest$Structure


# Two factor solution, promax rotation

EFAest <- fa(ECSICovData, nfactors = 2, rotate = "promax")

# Pattern coefficients are printed out by default
EFAest

# Structure coefficients differ from pattern coefficients
EFAest$Structure


##### Factor correlation (obtained from EFA) #####


# Four factor solution, promax rotation. Factor correlations are printed out
# by default

fa(ECSICovData, nfactors = 4, rotate = "promax")


#####  Techniques that require SEM software ##### 


##### Factor correlation (point estimate) ##### 

CFAest <- cfa(CFAmodel, sample.cov = ECSICovData, sample.nobs = N,
              std.lv = TRUE) # use alternative scaling

summary(CFAest)


##### Factor correlation (whether the confidence interval includes 1) ##### 

summary(CFAest, ci = TRUE)



##### Techniques using model fit indices: compared to ##### 
##### nested models with fewer factors                #####

CFAmodel2 <- "ACSIPERQ =~ ACSI1 + ACSI2 + ACSI3 + PERQ1 + PERQ2 + PERQ3
CUEX =~ CUEX1 + CUEX2 + CUEX3
PERV =~ PERV1 + PERV2"

CFAest2 <- cfa(CFAmodel2, sample.cov = ECSICovData, sample.nobs = N,
               std.lv = TRUE) # use alternative scaling

lavTestLRT(CFAest, CFAest2)

# Same analysis implemented with constrained correlations

CFAmodel3 <- "ACSI =~ ACSI1 + ACSI2 + ACSI3
CUEX =~ CUEX1 + CUEX2 + CUEX3
PERQ =~ PERQ1 + PERQ2 + PERQ3
PERV =~ PERV1 + PERV2

ACSI ~~ 1*PERQ
ACSI ~~ a*CUEX
PERQ ~~ a*CUEX
ACSI ~~ b*PERV
PERQ ~~ b*PERV"


CFAest3 <- cfa(CFAmodel3, sample.cov = ECSICovData, sample.nobs = N,
               std.lv = TRUE) # use alternative scaling

lavTestLRT(CFAest, CFAest3)



##### Techniques using model fit indices: comparison #####
##### against model with correlation fixed at 1      #####

CFAmodel4 <- "ACSI =~ ACSI1 + ACSI2 + ACSI3
CUEX =~ CUEX1 + CUEX2 + CUEX3
PERQ =~ PERQ1 + PERQ2 + PERQ3
PERV =~ PERV1 + PERV2

ACSI ~~ 1*PERQ"


CFAest4 <- cfa(CFAmodel4, sample.cov = ECSICovData, sample.nobs = N,
               std.lv = TRUE) # use alternative scaling

lavTestLRT(CFAest, CFAest4)

# All correlation pairs programmatically

Baseline <- CFAest

# Form the constraints
constraints <- paste(Baseline@ParTable$lhs[27:32], "~~1*",
                     Baseline@ParTable$rhs[27:32], sep="")

tests <- lapply(constraints,function(constraint){
   res <- lavTestLRT(Baseline,update(Baseline, add=constraint))
   rownames(res)[2] <- constraint # Add the model name
   res
})

do.call(rbind,c(tests[1],
                # Remove the baseline from all but first
                lapply(tests[-1], function(x){x[-1,]}))) 


##### Techniques using model fit indices: comparison #####
##### against model with correlation fixed at cutoff #####

CFAmodel5 <- "ACSI =~ ACSI1 + ACSI2 + ACSI3
CUEX =~ CUEX1 + CUEX2 + CUEX3
PERQ =~ PERQ1 + PERQ2 + PERQ3
PERV =~ PERV1 + PERV2

ACSI ~~ .96*PERQ"


CFAest5 <- cfa(CFAmodel5, sample.cov = ECSICovData, sample.nobs = N,
               std.lv = TRUE) # use alternative scaling

lavTestLRT(CFAest, CFAest5)


##### AVE: compared with the square of factor correlation #####

AVEs <- reliability(CFAest)["avevar",1:4]

comparisonMatrix <- inspect(CFAest, what = "cor.lv")
class(comparisonMatrix) <- "matrix"

diag(comparisonMatrix) <- AVEs

comparisonMatrix[upper.tri(comparisonMatrix)] <- 
  comparisonMatrix[upper.tri(comparisonMatrix)]^2

comparisonMatrix

##### AVE: compared with the square of scale score correlation #####

comparisonMatrix <- scaleScoreCorrelations

diag(comparisonMatrix) <- AVEs

comparisonMatrix[upper.tri(comparisonMatrix)] <- 
  comparisonMatrix[upper.tri(comparisonMatrix)]^2

comparisonMatrix


##### Structure coefficients (obtained from CFA) #####

inspect(CFAest, what = "cor.all")[1:11,12:15]

