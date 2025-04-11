## Import hand-labeled data
daily_LA_Female_drug <- read.csv("~/Desktop/Proposal/daily_LA_Female_drug.csv")
smp_size <- floor(0.67 * nrow(daily_LA_Female_drug))

## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(daily_LA_Female_drug)), size = smp_size)

train <- daily_LA_Female_drug[train_ind, ]
test <- daily_LA_Female_drug[-train_ind, ]

# check distribution of copycat/non-copycat states
copycat <- subset(train, Copycat.State == "C")
nocopy <- subset(train, Copycat.State == "NC")


x1 <- copycat$Frequency
fitp <- fitdist(x1, "pois")
fitn <- fitdist(x1, "nbinom")
cdfcomp(list(fitp, fitn), legendtext = c("Poisson", "negative binomial"), main = "Copycat")
# plot(fitp)
# plot(fitn)

x2 <- nocopy$Frequency
fitp <- fitdist(x2, "pois")
fitn <- fitdist(x2, "nbinom")
cdfcomp(list(fitp, fitn), legendtext = c("Poisson", "negative binomial"), main = "No Copycat")
# plot(fitp)
# plot(fitn)


# It actually shows mixture distribution
# d1 <- density(copycat$Frequency)
# d2 <- density(nocopy$Frequency)
# plot(d1, col = "red", main = "", xlab = "Counts of Suicide")
# lines(d2, col = "blue")
# legend(49, 0.01, c("Copycat", "No Copycat"), bty = "n", lty = c(1,1), 
#        col = c("red", "blue"))
# 
# 


sum(train[1:100,][3]=="C")/100
sum(train[1:100,][3]=="NC")/100

xtrain <- train$Frequency
strain <- train$Copycat.State # don't forget . between the name
Ntrain <- length(xtrain)
train <- list(x=xtrain, s = strain, N=Ntrain)
class(train) <- "hsmm.data"

xtest <- test$Frequency
stest <- test$Copycat.State # don't forget . between the name
Ntest <- length(xtest)
test <- list(x=xtest, s = stest, N=Ntest)
class(test) <- "hsmm.data"

p_copy <- summary(as.factor(train$s))[1]/(summary(as.factor(train$s))[1]+summary(as.factor(train$s))[2])
p_nocopy <- summary(as.factor(train$s))[2]/(summary(as.factor(train$s))[1]+summary(as.factor(train$s))[2])
trans_mat <- markovchainFit(strain)$estimate
trans_mat

j <- 2
initial <- c(p_copy, p_nocopy)
p <- matrix(c(0.05504587,0.04411135,0.9449541,0.9558887),nrow=j)
# mean of the two states in first 600 observations
b <- list(lambda = c(0.966,0.30))
model <- hmmspec(init=initial, trans=p, parms.emission=b,dens.emission=dpois.hsmm)
#plot(train, xlim = c(0,100))
h1 <- hmmfit(train,model,mstep=mstep.pois)
summary(h1)
ypred = predict(h1, test)
summary(as.factor(ypred$s))
ypred$s
ypred$x
accuracy1 <- ypred$s == as.numeric(stest)
accuracy1 <- as.data.frame(accuracy1)
accuracy1 <- as.factor(accuracy1$accuracy1)
accuracy <- 100*(sum(accuracy1 == "TRUE")/length(accuracy1))
accuracy
missclass <- 100*(sum(accuracy1 == "FALSE")/length(accuracy1))
missclass

#0: NC, 1: C
stest <- factor(stest, labels = c("0","1")) # substitute NC and C with 0 and 1 respectively
stest <- unfactor(stest)
pred <- ypred$s
pred <- factor(pred, labels = c("0","1")) # substituting 1 and 2 with 0 and 1 respectively.
pred <- unfactor(pred) 

tableMat <- table(predict=pred, truth=stest)
Specificity = tableMat[1,1]/(sum(tableMat[,1]))
Sensitivity = tableMat[2,2]/(sum(tableMat[,2]))
tableMat
Specificity
Sensitivity


# validate with synthetic data

## sythetic data validate

Synthetic_Data <- read_excel("~/Desktop/Proposal/Synthetic Data.xlsx")# will there be read.excel: sometimes
Synthetic_Data$states <- factor(Synthetic_Data$state,
                                labels = c("NC", "C"))
xtest <- Synthetic_Data$count
stest <- Synthetic_Data$states
Ntest <- length(xtest)
test <- list(x=xtest, N=Ntest)
class(test) <- "hsmm.data"

ypred = predict(h1, test)
# 1: NC, 2: C
#0: NC, 1: C
stest <- factor(stest, labels = c("0","1")) # substitute NC and C with 0 and 1 respectively
stest <- unfactor(stest)
pred <- ypred$s
pred <- factor(pred, labels = c("0","1")) # substituting 1 and 2 with 0 and 1 respectively.
pred <- unfactor(pred) 

tableMat <- table(predict=pred, truth=stest)
Specificity = tableMat[1,1]/(sum(tableMat[,1]))
Sensitivity = tableMat[2,2]/(sum(tableMat[,2]))
tableMat
Specificity
Sensitivity


#####################################
## Dr. Osgood Suggestion ############
#####################################

smp_size <- floor(0.67 * nrow(Synthetic_Data))

## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(Synthetic_Data)), size = smp_size)

train <- Synthetic_Data[train_ind, ]
test <- Synthetic_Data[-train_ind, ]


xtrain <- train$count
strain <- train$states
Ntrain <- length(xtrain)
train <- list(x=xtrain, N=Ntrain)
class(train) <- "hsmm.data"

xtest <- test$count
stest <- test$states
Ntest <- length(xtest)
test <- list(x=xtest, N=Ntest)
class(test) <- "hsmm.data"

p_copy <- summary(as.factor(train$states))[2]/(summary(as.factor(train$states))[1]+summary(as.factor(train$states))[2])
p_nocopy <- summary(as.factor(train$s))[1]/(summary(as.factor(train$s))[1]+summary(as.factor(train$s))[2])
trans_mat <- markovchainFit(strain)$estimate
trans_mat

j <- 2
initial <- c(0.5, 0.5)
p <- matrix(c(0.08427162, 0.91572849, 0.08604840, 0.9139516),nrow=j)
b <- list(lambda = c(2,0.01))
model <- hmmspec(init=initial, trans=p, parms.emission=b,dens.emission=dpois.hsmm)
#plot(train, xlim = c(0,100))
h1 <- hmmfit(train,model,mstep=mstep.pois)
summary(h1)
ypred = predict(h1, test)
summary(as.factor(ypred$s))
ypred$s
ypred$x
accuracy1 <- ypred$s == as.numeric(stest)
accuracy1 <- as.data.frame(accuracy1)
accuracy1 <- as.factor(accuracy1$accuracy1)
accuracy <- 100*(sum(accuracy1 == "TRUE")/length(accuracy1))
accuracy
missclass <- 100*(sum(accuracy1 == "FALSE")/length(accuracy1))
missclass

