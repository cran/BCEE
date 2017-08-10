ABCEE <- function(X, Y, U, omega, forX = NA, niter = 5000, nburn = 500, nthin = 10, maxmodelY = NA, OR = 20, family.X = "gaussian")
{
	n = length(Y); #Sample size
	n_cov = ncol(as.matrix(U)); #Number of covariate, potential confounders
	priorX = rep(0, n_cov);
	if(is.na(forX[1])) forX = 1:n_cov;
	priorX[forX] = 0.5;
	if(is.na(maxmodelY)) maxmodelY = min(niter + nburn, 2**n_cov);
	if(n_cov < 50) #bic.glm doesn't seem to work with > 50 covariates
	{
		model.X = bic.glm(y = X, x = U, glm.family = family.X, OR = OR, prior.param = priorX, maxCol = n_cov + 1);
		models.X = cbind(model.X$which, model.X$postprob);
		alpha_X = model.X$probne0/100;
	} else if((class(family.X) == "character" && family.X == "gaussian") |
		(class(family.X) == "family" & family.X$family == "gaussian" & family.X$link == "identity")) #Approximate posterior distribution of the exposure model, but 
	{
		resultsX = summary(regsubsets(y = X, x = U, nbest = 150, really.big = T, nvmax = n_cov))
		MLx = exp(-resultsX$bic/2 + min(resultsX$bic)/2);
		MLx = MLx/sum(MLx);
		models.X = cbind(resultsX$which[, -1], MLx);
		alpha_X = colSums(apply(models.X, 2, "*", MLx));
	} else
	{
		k_X = 1; #Index of pX_tested to be filled next 
		j = 1;
		resultsX = summary(regsubsets(y = X, x = U, nbest = 1, really.big = T, nvmax = n_cov)); #use linear regression for first screening of models
		alpha_X = resultsX$which[which.min(resultsX$bic), -1]; #Inital exposure model is the best fitting linear regression
	
		###Obtaining BIC for X
		if(sum(alpha_X) == 0) 
		{
			model_X = glm(X~1, family = family.X);
		} else
		{
			model_X = glm(X~U[,alpha_X == 1], family = family.X);
		}
		bic_x = BIC(model_X);

		models_X = matrix(0, nrow = ceiling(niter/nthin), ncol = n_cov); #objects that will contain results
		tested_models_X = matrix(-1, nrow = maxmodelY, ncol = 1); #objects that will contain information for models already tested
		pX_tested = matrix(NA, nrow = maxmodelY, ncol = 1); #objects that will contain information for models already tested

		for(i in 1:(nburn+niter))
		{
			alpha_X_0 = alpha_X; #Actual alpha_X
			alpha_X_1 = alpha_X; temp = sample(n_cov, 1); alpha_X_1[temp] = (alpha_X_1[temp] + 1)%%2; 
						#Candidate alpha_X
						#According to MC3, we select a model in the neighborhood of the candidate model,
						#that is, a model with one more variable or one variable fewer. Each model
						#have the same probability. Here, I randomly select one of the n_cov potential
						#confounder and either remove it or add it according to whether the variable
						#was already in the model or not, respectively.
			bic_x_0 = bic_x;
	
			char_alpha_X_1 = paste(alpha_X_1, collapse = "");
			tested = which(tested_models_X == char_alpha_X_1);
			if(length(tested) == 1) #The candidate model was already tested
			{
				bic_x_1 = pX_tested[tested];
			} else #The candidate model was never tested
			{
				if(k_X > maxmodelY)
				{
					print("Error: The number of models tested is larger than maxmodelY.");
					return();
				}
				###Obtaining elements to calculate B10 for X
				if(sum(alpha_X_1) == 0) 
				{
					model_X = glm(X~1, family = family.X);

				}
				else
				{
					model_X = glm(X~U[,alpha_X == 1], family = family.X);
				}
				bic_x_1 = BIC(model_X);
				tested_models_X[k_X] = char_alpha_X_1; #Recording info about candidate model
				pX_tested[k_X] = bic_x_1; #Recording info about candidate model
				k_X = k_X + 1; #Changing index
			}
		
			ratio = exp(-bic_x_1/2 + bic_x_0/2);
		
			if(sample(c(1,0), size = 1, prob = c(min(1, ratio), 1 - min(1, ratio))))
			{
				alpha_X = alpha_X_1;
				bic_x = bic_x_1;
			}
			if(i > nburn & (i - nburn)%%nthin == 1)
			{
				models_X[j, ] = alpha_X;
				j = j + 1;
			}
		} #end of loop
		alpha_X = colMeans(models_X);
		models.X = as.matrix(aggregate(rep(1, nrow(models_X)) ~., data = models_X,
			FUN = function(x){sum(x)/nrow(models_X)}));
		colnames(models.X) = NULL;
	}
	if(n_cov < 49) #bic.glm doesn't seem to work with > 50 covariates
	{
		alpha_Y = as.numeric(bic.glm(y = Y, x = cbind(X, U), glm.family = "gaussian", OR = 1.0000001,
			prior.param = c(1, rep(0.5, n_cov)), maxCol = n_cov + 2)$which[1, 2:(n_cov+1)]); #Initial outcome model
	} else  
	{
		resultsY = summary(regsubsets(y = Y, x = cbind(X, U), nbest = 1, really.big = T, nvmax = n_cov + 1, force.in = 1));
		alpha_Y = resultsY$which[which.min(resultsY$bic), -c(1,2)];
	}



	if(sum(alpha_Y) == 0)
	{
		model_Y0 = clm(Xarg = cbind(1, X), yarg = Y);
	}
	else
	{
		model_Y0 = clm(Xarg = as.matrix(cbind(1, X, U[ , alpha_Y == 1])), yarg = Y);
	}

	ll = 0.5 * (- n * (log(2 * pi) + 1 - log(n) + log(sum(model_Y0$res^2))));
	df.ll = sum(alpha_Y) + 3; # nb de param + 1 (3 = 1 + intercept + X)
	bic_y = -2 * ll + log(n) * df.ll;

	betas = numeric(ceiling(niter/nthin)); #objects that will contain results
	models_Y = matrix(0, nrow = ceiling(niter/nthin), ncol = n_cov); #objects that will contain results
	tested_models_Y = matrix(-1, nrow = maxmodelY, ncol = 1); #objects that will contain information for models already tested
	pY_tested = matrix(NA, nrow = maxmodelY, ncol = 1); #objects that will contain information for models already tested
	pY_tested1 = matrix(NA, nrow = maxmodelY, ncol = n_cov);
	pY_tested2 = matrix(NA, nrow = maxmodelY, ncol = n_cov);
	j = 1; #Index to be filled next
	k_Y = 1; #Index of pY_tested to be filled next
	tested_models_Y[k_Y] = paste0(alpha_Y, collapse = ""); #Recording info about candidate model
	sy = sd(Y);
	su = apply(U, 2, sd);

	if(sum(alpha_Y != 0))
	{
		a = coef(model_Y0)[-(1:2)];
		MSE <- mean((model_Y0$fitted.values - Y)^2);
		var_betaHat <- MSE * solve(t(cbind(1, X, U[, alpha_Y == 1])) %*% cbind(1, X, U[, alpha_Y == 1]))
		b = sqrt(diag(var_betaHat))[-(1:2)]
	}
	else
	{
		a = NA;
		b = NA;
	}
	pY_tested1[k_Y, 1:sum(alpha_Y)] = a;
	pY_tested2[k_Y, 1:sum(alpha_Y)] = b;
	pY_tested[k_Y] = bic_y;
	k_Y = k_Y + 1; #Changing index
	for(i in 1:(nburn+niter))
	{
		alpha_Y_0 = alpha_Y; #Actual alpha_Y
		alpha_Y_1 = alpha_Y; temp = sample(n_cov, 1); alpha_Y_1[temp] = (alpha_Y_1[temp] + 1)%%2; #Candidate alpha_Y
		bic_y_0 = bic_y;
		a_0 = a;
		b_0 = b;
		
		char_alpha_Y_1 = paste0(alpha_Y_1, collapse = "");
		tested = which(tested_models_Y == char_alpha_Y_1);
		if(length(tested)) #The candidate model was already tested
		{
			a_1 = pY_tested1[tested,];
			b_1 = pY_tested2[tested,];
			bic_y_1 = pY_tested[tested];
		}
		else #The candidate model was never tested
		{
			if(k_Y > maxmodelY)
			{
				print("Error: The number of models tested is larger than maxmodelY.");
				return();
			}
			if(sum(alpha_Y_1) == 0)
			{
				model_Y1 = clm(Xarg = cbind(1, X), yarg = Y);
				a_1 = NA;
				b_1 = NA;
			}
			else
			{
				model_Y1 = clm(Xarg = as.matrix(cbind(1, X, U[ , alpha_Y_1 == 1])), yarg = Y);
				a_1 = coef(model_Y1)[-(1:2)];
				MSE1 <- mean((model_Y1$fitted.values - Y)^2);
				var_betaHat1 <- MSE1 * solve(t(cbind(1, X, U[, alpha_Y_1 == 1])) %*% cbind(1, X, U[, alpha_Y_1 == 1]))
				b_1 = sqrt(diag(var_betaHat1))[-(1:2)];
			}
			pY_tested1[k_Y, 1:sum(alpha_Y_1) ] = a_1;
			pY_tested2[k_Y, 1:sum(alpha_Y_1) ] = b_1;
			ll1 = 0.5 * (- n * (log(2 * pi) + 1 - log(n) + log(sum(model_Y1$res^2))));
			df.ll1 = sum(alpha_Y_1) + 3;
			bic_y_1 = -2 * ll1 + log(n) * df.ll1;
			tested_models_Y[k_Y] = char_alpha_Y_1; #Recording info about candidate model
			pY_tested[k_Y] = bic_y_1; #Recording info about candidate model
			k_Y = k_Y + 1; #Changing index
		}
		change_Y = alpha_Y_1 - alpha_Y_0;
		px = alpha_X[change_Y != 0];
		add_Y = sum(change_Y); #Is it proposed to add or remove a variable
		if(add_Y == 1)
		{
			delta = omega*(rnorm(1, a_1[change_Y[alpha_Y_1 != 0] != 0], sd = b_1[change_Y[alpha_Y_1 != 0] != 0])*su[change_Y != 0]/sy)**2;
		}
		else
		{
			delta = omega*(rnorm(1, a_0[change_Y[alpha_Y_0 != 0] != 0],sd = b_0[change_Y[alpha_Y_0 != 0] != 0])*su[change_Y != 0]/sy)**2;
		}
		if(add_Y == 1)
		{
			ratio_P_alpha = (px*(delta/(1 + delta)) + (1 - px)*0.5)/(px*(1/(1 + delta)) + (1 - px)*0.5);
		}
		if(add_Y == -1)
		{
			ratio_P_alpha = (px*(1/(1 + delta)) + (1 - px)*0.5)/(px*(delta/(1 + delta)) + (1 - px)*0.5);
		}
		if(ratio_P_alpha == 0)
		{
			ratio = 0;
		}
		else if(ratio_P_alpha == Inf)
		{
			ratio = 1;
		}
		else
		{
			ratio = exp(-bic_y_1/2 + bic_y_0/2)*ratio_P_alpha;
		}
		if(sample(c(1,0), size = 1, prob = c(min(1, ratio), 1 - min(1, ratio))))
		{
			alpha_Y = alpha_Y_1;
			bic_y = bic_y_1;
			a = a_1;
			b = b_1;
		}
		if(i > nburn & (i - nburn)%%nthin == 1)
		{
			if(sum(alpha_Y) != 0)
			{
				model_b = clm(Xarg = cbind(1, X, U[,alpha_Y == 1]), yarg = Y);
				MSEb <- mean((model_b$fitted.values - Y)^2);
				var_betaHatb <- MSEb * solve(t(cbind(1, X, U[, alpha_Y == 1])) %*% cbind(1, X, U[, alpha_Y == 1]))
				b_b = sqrt(diag(var_betaHatb))[2];
			}
			else
			{
				model_b = clm(Xarg = cbind(1, X), yarg = Y);
				MSEb <- mean((model_b$fitted.values - Y)^2);
				var_betaHatb <- MSEb * solve(t(cbind(1, X)) %*% cbind(1, X))
				b_b = sqrt(diag(var_betaHatb))[2];
			}
			beta_alphaY = rnorm(1, mean = model_b$coefficients[2], sd = b_b);
			betas[j] = beta_alphaY;
			models_Y[j, ] = alpha_Y;
			j = j + 1;
		}
	}
return(list(betas = betas, models.X = models.X, models.Y = models_Y));
}
