# Rewriting the l2a method for an FLQuant
# Works nice and fast for iters up to about 500
# Then runs out of memory and goes mental

l2a_func <- function(object, model, stat="sum", weights=FLQuant(1, dimnames=dimnames(object)), ...){
	# constants
	cat("Converting lengths to ages ...\n")
	dnms <- dimnames(object)
	if(!all.equal(dnms, dimnames(weights))) stop("Weights must have the same dimensions as the data.")
	len <- as.numeric(dnms[[1]])
	mit <- niters(model)
	qit <- length(dnms[[6]]) 
	if(mit>1 & qit>1 & mit!=qit) stop("Can not operate with diferent iterations in each object.")

	# model
	mod <- FLModelSim(model=grInvMod(model), params=params(model), vcov=vcov(model), distr=distr(model))
	age <- floor(predict(mod, len=len))

	# aggregates all lengths above linf
	age <- apply(age, 2, function(x){
		x[which.max(x[!is.infinite(x)]):length(x)] <- max(x[!is.infinite(x)], na.rm=T)
		x
	})
	# slicing and aggregating

	# new flq
	dnms[[1]] <- as.character(sort(unique(c(age))))
	if(mit>qit) dnms[[6]] <- as.character(1:mit)
	flq <- FLQuant(NA, dimnames=dnms, quant="age")

	# loop :(
	# both have one iter
	if(mit>qit){
		object <- object[,,,,,rep(1, mit)]		
		dimnames(object)$iter <- 1:mit
		weights <- weights[,,,,,rep(1, mit)]
		dimnames(weights)$iter <- 1:mit
		}
	if(mit<qit){
		age <- age[,rep(1,qit)]
		dimnames(age)$iter <- 1:qit
		}
	

#browser()
##------------------------------
#iter <- 1
#ob_it <- iter(object,iter)
#dnms_iter <- c(list(age = age[,iter]), dimnames(ob_it)[2:6])
#dimnames(ob_it) <- dnms_iter



#-----------------------------

    # Using 'with' courtesy of Colin
    om <- melt(object)
    dimnames(age) <- list(len = len, iter = dimnames(object)$iter)
    # Make a single vector of ages from age, repeating each one by nyears, nseasons, nareas, nunits
    n <- prod(dim(flq)[c(2,3,4,5)])
    age_vec <- c(apply(age,2,rep,n))
    #length(age_vec) == dim(om)[1] # sanity check
    # Sort om so it is ordered by iter and lengths go 1,2,3,.... 
    om <- om[order(om$iter,om$area, om$season, om$unit, om$year,  om$len),]
    om <- cbind(om, age = age_vec)
    rm(age_vec)
    gc()

    if (stat == "sum"){
        # Now operate - depending on stat
        #flq_slice <- with(om, tapply(value, list(age, year, unit, season, area, iter), sum, na.rm = TRUE)) 
        flq_slice <- tapply(om$value, list(om$age, om$year, om$unit, om$season, om$area, om$iter), sum, na.rm = TRUE) 
        rm(om) # This can get too big and smash your system so delete it
        gc() # And call the garbage collector
        dimnames(flq_slice) <- dimnames(flq)
        flq <- FLQuant(flq_slice)
        flq[is.na(flq)] <- 0
    }

    if (stat == "mean"){
        weightm <- melt(weights)
        weightm <- weightm[order(weightm$iter,weightm$area, weightm$season, weightm$unit, weightm$year,  weightm$len),]
        om <- cbind(om, weight = weightm$value)
        scaled_weight <- with(om, tapply(weight, list(age, year, unit, season, area, iter), function(x){ x / sum(x)}))
        rm(weightm, weights)
        gc()
        om <- cbind(om, scaled_weight=unlist(scaled_weight))
        rm(scaled_weight)
        gc()
        flq_slice <- with(om, tapply(value*scaled_weight, list(age, year, unit, season, area, iter), sum, na.rm = TRUE)) 
        rm(om) # This can get too big and smash your system so delete it
        gc() # And call the garbage collector
        dimnames(flq_slice) <- dimnames(flq)
        flq <- FLQuant(flq_slice)
    }

	units(flq) <- units(object)
	return(flq)
}

#----------------------------------------------------------------------------------------------------
