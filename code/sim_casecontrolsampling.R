##  OVERVIEW
##  A population has a binary exposure E, a binary covariate FEM, and
##    a binary outcome D. 
##  This script takes multiple case-control samples from the population,
##    computes various summary statistics on the sample and compares to the population.

library(DescTools)
library(ggplot2)
library(dplyr)
library(gridExtra)
setwd("P://Teach//536-2024//Lecture1//")

#######################################
##          Functions                ##
#######################################
## function to turn p into logit(p)
logit <- function(p){log(p/(1-p))}
## function to transform a log odds back to a probability
expit <- function(x){exp(x)/(1 + exp(x))}

## take a case-control sample from a population 'pop'
## with a 0/1 variable D for controls/cases
## default is 100 cases and 100 controls
takeCCsample <- function(pop, ncase=100, ncntl=100){
  allcntls <- subset(pop, D==0)
  ncntls <- nrow(allcntls)
  cntls.sampled <- allcntls[sample(ncntls, size=ncase), ]
  allcases <- subset(pop, D==1)
  ncases <- nrow(allcases)
  cases.sampled <- allcases[sample(ncases, size=ncase), ]
  mysample <- rbind(cntls.sampled, cases.sampled)
  return(mysample)
}

## calculate some summaries of E, FEM, D data 
my.summary <-function(mydata){
  proportion.exposed <- mean(mydata$E)
  proportion.female <- mean(mydata$FEM)
  subset.exposed <- subset(mydata, E==1)
  subset.unexpos <- subset(mydata, E==0)
  subset.Diseased <- subset(mydata, D==1)
  subset.NotDis   <- subset(mydata, D==0)
  #Probability of D among Exposed
  pDgExp <- mean(subset.exposed$D)
  #Probability of D among UnExposed
  pDgUex <- mean(subset.unexpos$D)
  #Probability of Exposure among Diseased
  pEgD    <- mean(subset.Diseased$E)
  #Probability of Exposure among not-Diseased
  pEgDbar <- mean(subset.NotDis$E)
  #Relative Risk
  rr.E <- pDgExp / pDgUex
  #Convert probabilities to odds and calculate OR
  oddsDgExp <- pDgExp / (1-pDgExp)
  oddsDgUex <- pDgUex / (1-pDgUex)
  or.E <- oddsDgExp / oddsDgUex
  #logistic regression of disease on Exposure and FEM
  ms <- glm(D ~ E + FEM, family="binomial", data=mydata)
  beta_E   <- ms$coefficients[2]
  beta_FEM <- ms$coefficients[3]
  return(c(pE=proportion.exposed, pFEM=proportion.female, pDgE=pDgExp, pDgU=pDgUex, pEgD=pEgD, pEgNotD=pEgDbar, RR=rr.E, OR=or.E, betaE=beta_E, betaFEM=beta_FEM))
}

## Take a single case-control sample and calculate summary
do.one.sample.and.summarize <- function(population, ncase=100, ncntl=100){
  mysample <- takeCCsample(population, ncase=ncase, ncntl=ncntl)
  summary.statistics <- my.summary(mysample)
}

######################################
## Enter population Data as a table ##
## and turn into a data frame       ##
######################################

population.table <- data.frame(
  FEM=c(rep(0, 4), rep(1, 4)), 
  E=rep(c(0,0,1,1), 2),  
  D=rep(c(0,1), 4), 
  n.sub=c(47070, 310, 2428, 32, 47190, 444, 2483, 43))
mypop <- Untable(population.table, freq="n.sub")
summary.pop <- my.summary(mypop)

#Quick Look at frequency of D in subgroups of the population
population.outcomes <- c(summary.pop[4], 
                         summary.pop[3], 
                         mean(subset(mypop, FEM==0)$D), 
                         mean(subset(mypop, FEM==1)$D))
barplot(height=population.outcomes, names=c(expression(P(D*"|"* bar(E))), expression(P(D*"|"* E)), expression(P(D*"|"* MALE)), expression(P(D*"|"* FEM))), main="Frequency of D in Subgroups", space=c(0.1,0.1,0.3, 0.1))
# names = c("P(D|not E)", "P(D|E)", "P(D|male)", "P(D|female)")
# ggplot(data = data.frame(x = names, y = population.outcomes), 
       # mapping = aes(x, y)) +
  # geom_col() +
  # theme_bw() +
  # labs(title = "Frequency of D in Subgroups", x = "", y = "")



###########################################################################
##  Take many Case-Control Samples, calculate summaries of each sample   ##
##  and store in an object called "results"                              ##
###########################################################################
set.seed(1888)
results <- replicate(200, do.one.sample.and.summarize(mypop, ncase=100, ncntl=100))
results <- t(results)

############################################################################
##  Summary of Case-Control Samples Vs. Population - DISPLAY 1
############################################################################
par(mfrow = c(3, 2))
hist(results[, 1], main="P(E)", xlab="probabiliy of being exposed")
abline(v = summary.pop[1], col="darkblue", lwd=4)
hist(results[, 2], main="P(female)", xlab="probability of being female")
abline(v = summary.pop[2], col="darkblue", lwd=4)

hist(results[, 3], main="P(D|E)", xlab="Prob of D among exposed", xlim=c(0, max(results[, 3]+0.03)))
abline(v = summary.pop[3], col="darkblue", lwd=4)
hist(results[, 4], main=expression(P(D*"|"* bar(E))), xlab="Prob of D among unexposed", xlim=c(0, max(results[, 4]+0.03)))
abline(v = summary.pop[4], col="darkblue", lwd=4)

hist(results[, 5], main="P(E|D)", xlab="probabiliy of E among Diseased")
abline(v = summary.pop[5], col="darkblue", lwd=4)
hist(results[, 6], main=expression(P(E*"|"* bar(D))), xlab="probability of E among not D")
abline(v = summary.pop[6], col="darkblue", lwd=4)

# this way is way worse!
results <- data.frame(results)
p1 <- ggplot(results, aes(pE)) + 
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept = summary.pop[1], col = "darkblue", lwd = 1) +
  labs(title = "P(E)", x = "probability of being exposed") + theme_bw()
        
p2 <- ggplot(results, aes(pFEM)) + 
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=summary.pop[2], col = "darkblue", lwd = 1) +
  labs(title = "P(female)", x = "probability of being female") + theme_bw()

p3 <- ggplot(results, aes(pDgE)) + 
  geom_histogram(binwidth = 0.1) + 
  geom_vline(xintercept=summary.pop[3], col = "darkblue", lwd = 1) +
  labs(title = "P(D|E)", x = "Prob of D among exposed") + theme_bw()

p4 <- ggplot(results, aes(pDgU)) + 
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=summary.pop[4], col = "darkblue", lwd = 1) +
  labs(title = "P(D|not E)", x = "Prob of D among unexposed") + theme_classic()

p5 <- ggplot(results, aes(pEgD)) + 
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=summary.pop[5], col = "darkblue", lwd = 1) +
  labs(title = "P(E|D)", x = "probabiliy of E among Diseased") + theme_bw()

p6 <- ggplot(results, aes(pEgNotD)) + 
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=summary.pop[6], col = "darkblue", lwd = 1) +
  labs(title = "P(E|not D)", x = "probability of E among not D") + theme_bw()

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)



## Each sample contained 100 randomly sampled cases and 100 randomly sampled controls
## Q1  Does the frequency of E in the sample  estimate the frequency of E in the population?  
##     In other words, do case-control samples provide an unbiased estimate of the frequency of E in the population?

## Response: The distribution of P(E) from the set of balanced case-control 
## samples slightly overestimates the true prevalence of exposure.

## Q2  Does the frequency of females in the sample estimate the frequency of females in the population?  
##     If not, can you reason why they are different?

## Response: The distribution of proportions from the set of balanced 
## case-control samples slightly overestimates the true proportion of female individuals.

## Q3 Does the frequency of D among Exposed in a case-control study estimate the frequency in the population?
##    Among Unexposed?

## Response: The distribution of frequencies from the set of balanced 
## case-control sample overestimates the true frequency of disease within 
## exposed individuals. The true frequency within unexposed individuals is 
## severely underestimated by the case-control samples.

## Q4 Does the frequency of E among diseased in a case-control study estimate the frequency in the population?
##    Does the frequency of E among non-diseased in a case-control study estimate the frequency in the population?

## Response: The 

## Q5 Compare your answers for Q3 and Q4.  What explains the difference?



############################################################################
##  Summary of Case-Control Samples Vs. Population - DISPLAY 2
############################################################################

par(mfrow = c(3, 2))
hist(results[, 7], main="RR, exposed vs. unexposed", xlab="RR", sub="")
abline(v = summary.pop[7], col="darkblue", lwd=4)
hist(results[, 8], main="OR, exposed vs. unexposed", xlab="OR", sub="")
abline(v = summary.pop[8], col="darkblue", lwd=4)


hist(results[, 9], main="regr coeff for E", xlab="log(OR)")
abline(v = summary.pop[9], col="darkblue", lwd=4)
hist(results[, 10], main="regr coeff for FEM", xlab="log(OR)")
abline(v = summary.pop[10], col="darkblue", lwd=4)


## Each sample contained 100 randomly sampled cases and 100 randomly sampled controls
## Q6  The first plot shows the relative risk of D comparing exposed and unexposed.
##     Does this quantity calculated in case-control samples estimate the RR in the population?

## Q7  The second plot shows the Odds Ratio instead of the relative risk.
##     Does this quantity calculated in case-control samples estimate the OR in the population?

## Q8  'reg logOR E' is the coefficient of E in a logistic regression of D on E and FEM
##      Does the quantity in case-control samples estimate the quantity in the population?
##      Consider also 're logOR FEM', which is the coefficient for FEM in the same model


