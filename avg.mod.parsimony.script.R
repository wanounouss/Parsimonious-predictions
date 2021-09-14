##
## do we need a good nechanistic understanding of the generating model to make accurate predictions?
## what model to use to make predictions? averaged model? full model? reduced model?
## do inferences about which variables are important differ when making them based on full models compared to averaged models?
## is it a problem to include confonding variables (whose effects taper) in models for inferences?
## does the model's complexcoty affect inferences about the generating phenomenun (i.e. which variable to discard?)
##



# install.packages("MuMIn")
library(MASS)
library(MuMIn)

#	Dataset simulation function - 5 genuine + 2 spurious explanatory variables

dat.sim.5var <- function(SampleSize=SampleSize, r1=r1, r2=r2, r3=r3, r6=r6, r5=r5, r1_2=r1_2) {

	# correlation matrix 
	M = matrix(c( 1, r1  , r2   , r3 ,  r6 , r5, 
				  r1,  1   , r1_2 , 0  ,  0  , 0 , 
				  r2,  r1_2,  1 , 0  ,  0  , 0 ,
				  r3,  0   ,  0   , 1,  0  , 0 ,
				  r6,  0   ,  0   , 0  , 1 , 0 ,
				  r5,  0   ,  0   , 0  ,  0  , 1
				 ), nrow=6, ncol=6)

		
	dat <- mvrnorm(n=SampleSize, mu=rep(0,6), Sigma=M)
	# dat <- mvrnorm(n=SampleSize, mu=c(30, rep(0, 5)), Sigma=M)
	dat <- as.data.frame(dat); names(dat) <- c("Y", "x1", "x2", "x3", "x6", "x5")
		
	dat$x4 <- rnorm(n=SampleSize)	# spurious variable
	dat$x7 <- rnorm(n=SampleSize)	# spurious variable
	
	return(dat)
}	
r.1 <- 0.6
r.2 <- 0.2
r.3 <- 0.05
r.6 <- 0.5
r.5 <- -0.4
r1.2 <- 0


##
## cross validation
##

na.func <- function(r){
	if(is.na(r) == FALSE){
		return(r)
	} else{
		return(0)
	}
}
thresh.Pval <- 0.05
thresh.SW <- 0.5
sample.size.vec <- c(20, 30, 50, 80, 100, 200, 500, 1000)
mse.first.sample.size <- numeric(length(sample.size.vec))
mse.tot.sample.size <- numeric(length(sample.size.vec))
mse.full.sample.size <- numeric(length(sample.size.vec))
mse.nat.sample.size <- numeric(length(sample.size.vec))
mse.reduced.full.sample.size <- numeric(length(sample.size.vec))
mse.reduced.nat.sample.size <- numeric(length(sample.size.vec))
mse.reduced.SW.full.sample.size <- numeric(length(sample.size.vec))
mse.reduced.SW.nat.sample.size <- numeric(length(sample.size.vec))
count <- 1
for(sample.size in sample.size.vec){ # sample.size <- 50
  print(sample.size)
	iteration <- 100
	mse.first <- numeric(iteration)
	mse.tot <- numeric(iteration)
	mse.full <- numeric(iteration)
	mse.nat <- numeric(iteration)
	mse.reduced.full <- numeric(iteration)
	mse.reduced.nat <- numeric(iteration)
	mse.reduced.SW.full <- numeric(iteration)
	mse.reduced.SW.nat <- numeric(iteration)
	for(i in 1:iteration){
		dtot <- dat.sim.5var(sample.size, r.1, r.2, r.3, r.6, r.5, r1.2)
		options(na.action = "na.fail") 
		d.train <- dtot[1:(sample.size/2),]
		d.test <- dtot[(sample.size/2+1):sample.size,]
		m <- lm(Y ~ x1 + x2 + x3 + x4, data = d.train)
		ms1 <- dredge(m, extra = "R^2", rank = "AIC")
		# ms1 <- dredge(m, extra = "R^2", rank = "AIC", m.lim = c(1, 1)) ## limits the number of parameter to be considered in models: e.g. a maximum of 2
		mod.avg <- model.avg(ms1)
		sum.mod.avg <- summary(mod.avg)
		imp <- importance(ms1)

		## first ranked model
		int.first <- ms1[1,1]
		r.1.first <- na.func(ms1$x1[1])
		r.2.first <- na.func(ms1$x2[1])
		r.3.first <- na.func(ms1$x3[1])
		r.4.first <- na.func(ms1$x4[1])
		
		## full model
		totm <- coefficients(m)
		int.tot <- totm[1]
		r.1.tot <- totm[2]
		r.2.tot <- totm[3]
		r.3.tot <- totm[4]
		r.4.tot <- totm[5]
		
		## averaged model
		full <- sum.mod.avg$coefmat.full
		nat <- sum.mod.avg$coefmat.subset 

		int.full <- full["(Intercept)", "Estimate"]
		r.1.full <- full["x1", "Estimate"]
		r.2.full <- full["x2", "Estimate"]
		r.3.full <- full["x3", "Estimate"]
		r.4.full <- full["x4", "Estimate"]

		int.nat <- nat["(Intercept)", "Estimate"]
		r.1.nat <- nat["x1", "Estimate"]
		r.2.nat <- nat["x2", "Estimate"]
		r.3.nat <- nat["x3", "Estimate"]
		r.4.nat <- nat["x4", "Estimate"]
		
		## reduced averaged model based on p values
		reduced.full <- full
		reduced.full[, "Estimate"][reduced.full[, "Pr(>|z|)"] >= thresh.Pval] <- 0
		reduced.nat <- nat
		reduced.nat[, "Estimate"][reduced.nat[, "Pr(>|z|)"] >= thresh.Pval] <- 0
		
		int.reduced.full <- reduced.full["(Intercept)", "Estimate"]
		r.1.reduced.full <- reduced.full["x1", "Estimate"]
		r.2.reduced.full <- reduced.full["x2", "Estimate"]
		r.3.reduced.full <- reduced.full["x3", "Estimate"]
		r.4.reduced.full <- reduced.full["x4", "Estimate"]

		int.reduced.nat <- reduced.nat["(Intercept)", "Estimate"]
		r.1.reduced.nat <- reduced.nat["x1", "Estimate"]
		r.2.reduced.nat <- reduced.nat["x2", "Estimate"]
		r.3.reduced.nat <- reduced.nat["x3", "Estimate"]
		r.4.reduced.nat <- reduced.nat["x4", "Estimate"]
		
		## reduced averaged model based on SW
		imp <- as.vector(sum.mod.avg$importance)
		reduced.SW.full <- full
		reduced.SW.full[, "Estimate"][imp <= thresh.SW] <- 0
		reduced.SW.nat <- nat
		reduced.SW.nat[, "Estimate"][imp <= thresh.SW] <- 0
		
		int.reduced.SW.full <- reduced.SW.full["(Intercept)", "Estimate"]
		r.1.reduced.SW.full <- reduced.SW.full["x1", "Estimate"]
		r.2.reduced.SW.full <- reduced.SW.full["x2", "Estimate"]
		r.3.reduced.SW.full <- reduced.SW.full["x3", "Estimate"]
		r.4.reduced.SW.full <- reduced.SW.full["x4", "Estimate"]

		int.reduced.SW.nat <- reduced.SW.nat["(Intercept)", "Estimate"]
		r.1.reduced.SW.nat <- reduced.SW.nat["x1", "Estimate"]
		r.2.reduced.SW.nat <- reduced.SW.nat["x2", "Estimate"]
		r.3.reduced.SW.nat <- reduced.SW.nat["x3", "Estimate"]
		r.4.reduced.SW.nat <- reduced.SW.nat["x4", "Estimate"]
		
		##
		## predictions
		##

		pred.first <- int.first + r.1.first*d.test$x1 + r.2.first*d.test$x2 + r.3.first*d.test$x3 + r.4.first*d.test$x4
		mse.first[i] <- mean((d.test$Y - pred.first)^2)
	
		pred.tot <- int.tot + r.1.tot*d.test$x1 + r.2.tot*d.test$x2 + r.3.tot*d.test$x3 + r.4.tot*d.test$x4
		mse.tot[i] <- mean((d.test$Y - pred.tot)^2)
		
		pred.full <- int.full + r.1.full*d.test$x1 + r.2.full*d.test$x2 + r.3.full*d.test$x3 + r.4.full*d.test$x4
		mse.full[i] <- mean((d.test$Y - pred.full)^2)

		pred.nat <- int.nat + r.1.nat*d.test$x1 + r.2.nat*d.test$x2 + r.3.nat*d.test$x3 + r.4.nat*d.test$x4
		mse.nat[i] <- mean((d.test$Y - pred.nat)^2)
		
		pred.reduced.full <- int.reduced.full + r.1.reduced.full*d.test$x1 + r.2.reduced.full*d.test$x2 + r.3.reduced.full*d.test$x3 + r.4.reduced.full*d.test$x4
		mse.reduced.full[i] <- mean((d.test$Y - pred.reduced.full)^2)

		pred.reduced.nat <- int.reduced.nat + r.1.reduced.nat*d.test$x1 + r.2.reduced.nat*d.test$x2 + r.3.reduced.nat*d.test$x3 + r.4.reduced.nat*d.test$x4
		mse.reduced.nat[i] <- mean((d.test$Y - pred.reduced.nat)^2)
		
		pred.reduced.SW.full <- int.reduced.SW.full + r.1.reduced.SW.full*d.test$x1 + r.2.reduced.SW.full*d.test$x2 + r.3.reduced.SW.full*d.test$x3 + r.4.reduced.SW.full*d.test$x4
		mse.reduced.SW.full[i] <- mean((d.test$Y - pred.reduced.SW.full)^2)

		pred.reduced.SW.nat <- int.reduced.SW.nat + r.1.reduced.SW.nat*d.test$x1 + r.2.reduced.SW.nat*d.test$x2 + r.3.reduced.SW.nat*d.test$x3 + r.4.reduced.SW.nat*d.test$x4
		mse.reduced.SW.nat[i] <- mean((d.test$Y - pred.reduced.SW.nat)^2)
		
	}
	mse.first.sample.size[count] <- mean(mse.first)
	mse.tot.sample.size[count] <- mean(mse.tot)
	mse.full.sample.size[count] <- mean(mse.full)
	mse.nat.sample.size[count] <- mean(mse.nat)
	mse.reduced.full.sample.size[count] <- mean(mse.reduced.full)
	mse.reduced.nat.sample.size[count] <- mean(mse.reduced.nat)
	mse.reduced.SW.full.sample.size[count] <- mean(mse.reduced.SW.full)
	mse.reduced.SW.nat.sample.size[count] <- mean(mse.reduced.SW.nat)
	count <- count + 1
}

plot(NULL, xlim = c(1, length(sample.size.vec)), ylim = c(0.6, max(c(mse.first.sample.size, mse.full.sample.size, mse.nat.sample.size, mse.tot.sample.size))),
	xlab = "sample size", ylab = "mse", pch = 19, xaxt = "n", las = 1, cex.lab = 1.5, cex.axis = 1.3)
points(mse.first.sample.size ~ as.factor(sample.size.vec), pch = 19)
points(mse.tot.sample.size ~ as.factor(sample.size.vec), pch = 17)
points(mse.full.sample.size ~ as.factor(sample.size.vec), pch = 15, col = "red")
points(mse.nat.sample.size ~ as.factor(sample.size.vec), col = "red")
points(mse.reduced.full.sample.size ~ as.factor(sample.size.vec), pch = 15, col = "blue")
points(mse.reduced.nat.sample.size ~ as.factor(sample.size.vec), col = "blue")
points(mse.reduced.SW.full.sample.size ~ as.factor(sample.size.vec), pch = 15, col = "green")
points(mse.reduced.SW.nat.sample.size ~ as.factor(sample.size.vec), col = "green")
axis(side = 1, at = c(1:8), labels = levels(as.factor(sample.size.vec/2)))
text(6, 1.25, "first ranked model")
text(6, 1.2, "full model")
text(6, 1.15, "full averaging")
text(6, 1.1, "natural averaging")
text(6, 1.05, "reduced full averaging P value")
text(6, 1, "reduced natural averaging P value")
text(6, 0.95, "reduced full averaging SW")
text(6, 0.9, "reduced natural averaging SW")
points(4.2, 1.25, pch = 19)
points(4.2, 1.2, pch = 17)
points(4.2, 1.15, pch = 15, col = "red")
points(4.2, 1.1, pch = 1, col = "red")
points(4.2, 1.05, pch = 15, col = "blue")
points(4.2, 1, pch = 1, col = "blue")
points(4.2, 0.95, pch = 15, col = "green")
points(4.2, 0.9, pch = 1, col = "green")










full <- summary(mod.avg)$coefmat.full
nat <- summary(mod.avg)$coefmat.subset 
iteration <- 1
AICfull <- numeric(length(iteration))
AICnat <- numeric(length(iteration))
for(i in 1:iteration){
	r.1.full <- rnorm(1, full["x1", "Estimate"], full["x1", "Std. Error"])
	r.2.full <- rnorm(1, full["x2", "Estimate"], full["x2", "Std. Error"])
	r.3.full <- rnorm(1, full["x3", "Estimate"], full["x3", "Std. Error"])
	r.5.full <- rnorm(1, full["x4", "Estimate"], full["x4", "Std. Error"])
	r.6 <- 0.05
	r1.2 <- 0

	dfull <- dat.sim.5var(sample.size, r.1.full, r.2.full, r.3.full, r.6, r.5.full, r1.2)
	options(na.action = "na.fail") 
	m.full <- lm(Y ~ x1 + x2 + x3 + x5, data = dfull)
	AICfull[i] <- AIC(m.full)

	r.1.nat <- rnorm(1, nat["x1", "Estimate"], nat["x1", "Std. Error"])
	r.2.nat <- rnorm(1, nat["x2", "Estimate"], nat["x2", "Std. Error"])
	r.3.nat <- rnorm(1, nat["x3", "Estimate"], nat["x3", "Std. Error"])
	r.5.nat <- rnorm(1, nat["x4", "Estimate"], nat["x4", "Std. Error"])
	r.6 <- 0.05
	r1.2 <- 0
	dnat <- dat.sim.5var(sample.size, r.1.nat, r.2.nat, r.3.nat, r.6, r.5.nat, r1.2)
	options(na.action = "na.fail") 
	m.nat <- lm(Y ~ x1 + x2 + x3 + x5, data = dnat)
	AICnat[i] <- AIC(m.nat)
}
ms1$AIC[1]
mean(AICfull) - ms1$AIC[1]
mean(AICnat) - ms1$AIC[1]
 












