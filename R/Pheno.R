CheckArguments <- function(x, ...) UseMethod("CheckArguments")

CheckArguments.Pheno <- function(object){
	
	data <- object$data 
	x <- object$x
	y <- object$y
	label <- object$label 
	defaultLabel <- object$defaultLabel
	block <- object$block
	ncluster <- object$ncluster
	orderby <- object$orderby
	method <- object$method
	step <- object$step
	width <- object$width
	nfeatures <- object$nfeatures
	
	# check data argument
	if(is.numeric(nrow(data)) != TRUE)
		stop("invalid 'data'")
	if(nrow(na.omit(data)) == 0)
		stop("invalid 'data' after removing NA vaule")
	
	# check x, y arguments
	if((is.null(x) | is.null(y)) == TRUE)
		stop("invalid 'x' or 'y'")
	
	# check names, replicated names are allowed by using [[]]
	varNames <- colnames(data)
	inputNames <- c(x, y, label, block, orderby)
	if(all(inputNames %in% varNames) != TRUE)
		stop("invalid argument names")
	
	# check label and defaultLabel arguments
	if(!is.null(label) & is.null(defaultLabel))
		stop("invalid 'defaultLabel' but 'label' exists")
	if(!is.null(defaultLabel) & is.null(label))
		stop("invalid 'label' but 'defaultLabel' exists")
	if(!is.null(label) & length(label) >= 2)
		stop("invalid 'label', at most 1 dimension")
	if(!is.null(label)){
		if(length(unique(data[[label]])) < 2)
			stop("invalid 'label', which at least has two different components")
		if(!is.null(defaultLabel) & length(defaultLabel) >= length(unique(data[[label]])))
			stop("invalid 'defaultLabel', it must have less components than 'label'")	
		if(!is.null(defaultLabel) & !(defaultLabel %in% unique(data[[label]])))
			stop("invalid 'defaultLabel', which should be any component of 'label'")
	}	
	
	# check ncluster argument
	if(is.null(label) & is.null(ncluster))
		stop("invalid 'ncluster', the number of clusters")
	
	# check pcNumber argument
#	if(iis.null(pcNumber))
#		stop("invalid 'pcNumber', the number of principle components")
#	if(!(pcNumber == round(pcNumber)) | (pcNumber < 0))
#		stop("invalid 'pcNumber', it must be a positive integer")
		
	# check block argument
	if(!is.null(block) & length(block) > 2)
		stop("invalid 'block', at most 2 dimensions")
	
	# check orderby argument
	if(!is.null(orderby) & length(orderby) >= 2)
		stop("invalid 'orderby', at most 1 dimension")
	
	# check method argument	
	if(is.null(method) | length(method) != 1)
		stop("invalid 'method', at most 1 dimension")
	if(!(method %in% c("SVM", "RF", "KMEANS")))
		stop("invalid 'method', only 'SVM', 'RF', 'KMEANS' are available")
	if(!is.null(label) & (method %in% c("KMEANS")))
		stop("invalid 'method', must be 'SVM' or 'RF' due to existing 'label'")
	if(is.null(label) & (method %in% c("SVM", "RF")))
		stop("invalid 'method', must be 'KMEANS' due to missing 'label'")
	
	# check step, width, nfeatures arguments
	if(!(step == round(step)) | (step < 0))
		stop("invalid 'step', it must be a positive integer")
	if(!(width == round(width)) | (width < 0))
		stop("invalid 'width', it must be a positive integer")
	if(!(nfeatures == round(nfeatures)) | (nfeatures < 0))
		stop("invalid 'nfeatures', it must be a positive integer")
}

CheckArguments.predictPheno <- function(object){
	
	data <- object$data 
	x <- object$x
	y <- object$y
	block <- object$block
	orderby <- object$orderby
	step <- object$step
	width <- object$width
	
	# check data argument
	if(is.numeric(nrow(data)) != TRUE)
		stop("invalid 'data'")
	if(nrow(na.omit(data)) == 0)
		stop("invalid 'data' after removing NA vaule")
	
	# check x, y arguments
	if((is.null(x) | is.null(y)) == TRUE)
		stop("invalid 'x' or 'y'")
	
	# check names, replicated names are allowed by using [[]]
	varNames <- colnames(data)
	inputNames <- c(x, y, block, orderby)
	if(all(inputNames %in% varNames) != TRUE)
		stop("invalid argument names")
	
	# check block argument
	if(!is.null(block) & length(block) > 2)
		stop("invalid 'block', at most 2 dimensions")
	
	# check orderby argument
	if(!is.null(orderby) & length(orderby) >= 2)
		stop("invalid 'orderby', at most 1 dimension")
	
	# check step, width, nfeatures arguments
	if(!(step == round(step)) | (step < 0))
		stop("invalid 'step', it must be a positive integer")
	if(!(width == round(width)) | (width < 0))
		stop("invalid 'width', it must be a positive integer")
	
}

TransferData <- function(data, x, y, label, block, orderby){
	
	# to data frame
	WorkingData <- as.data.frame(data)
	WorkingData <- na.omit(WorkingData)
	
	# block argument
	if(is.null(block)){
		WorkingData$blockTemp <- seq(nrow(WorkingData))
		block <- "blockTemp"
	} else{
		if(length(block) == 2){
			WorkingData$blockTemp <- paste(WorkingData[[block[1]]], 
				WorkingData[[block[2]]], sep = "")
			block <- "blockTemp"
		}
	}
	
	# orderby argument
	if(is.null(orderby)){
		warning("missing 'orderby', data is ordered by default")
		WorkingData$orderbyTemp <- seq(nrow(WorkingData))
		orderby <- "orderbyTemp"
	}
	
	# select subset data
	if(!is.null(label)){
		WorkingDataSub <- subset(WorkingData, select = c(
			x, y, label, block, orderby))
		WorkingDataSub <- WorkingDataSub[order(WorkingDataSub[[orderby]], 
			decreasing = FALSE), ]
	} else{
		WorkingDataSub <- subset(WorkingData, select = c(
			x, y, block, orderby))
		WorkingDataSub <- WorkingDataSub[order(WorkingDataSub[[orderby]], 
			decreasing = FALSE), ]
	}
	
	return(list(data = WorkingDataSub, block = block, orderby = orderby))
}

SplitData <- function(data, block, orderby, 
	testBlockProp, blockOrderedLabels, blockOrderedNames){
	
	blockDataFrame <- data.frame(Name = blockOrderedNames, 
		Label = blockOrderedLabels)
	
	splitDatabyLabels <- split(blockDataFrame, blockDataFrame$Label)
	countSplitDatabyLabels <- sapply(splitDatabyLabels, nrow, simplify = TRUE)
	countNamesTemp <- names(countSplitDatabyLabels)
	countLabelsTemp <- as.vector(countSplitDatabyLabels)
	testSizes <- ceiling(countLabelsTemp * testBlockProp)
	
	FUNSAMPLE <- function(x, y){
		return(sample(x, size = y, replace = FALSE,
			prob = rep(1 / x, x)))
	}
	testSample <- mapply(FUNSAMPLE, x = countLabelsTemp, y = testSizes, 
		SIMPLIFY = FALSE)
	names(testSample) <- countNamesTemp
	testNames <- c()
	for(i in seq(length(countNamesTemp))){
		testNames <- c(testNames, 
			splitDatabyLabels[[i]][testSample[[i]], ]$Name)
	}
	trainNames <- seq(nrow(blockDataFrame))[- testNames]
	testNames <- blockDataFrame$Name[mixedsort(testNames)]
	trainNames <- blockDataFrame$Name[mixedsort(trainNames)]
	
	testData <- data[which(data[[block]] %in% testNames), ]
	testData <- testData[order(testData[[orderby]]), ]
	trainData <- data[which(data[[block]] %in% trainNames), ]
	trainData <- trainData[order(trainData[[orderby]]), ]
	
	returnData <- list(test = testData, train = trainData, 
		testName = testNames)	
	return(list(data = returnData))
}


BAYNIGFUN <- function(x, y, step, width) {
	# compute Bayesian NIG for any input without spliting
	R2FUN <- function(y, Est.y) return(1 - sum((y - Est.y) ^ 2) / sum((y - mean(y)) ^ 2))
	BETAFUN <- function(x, y) return(lm(y ~ x)$coefficients)
	BETAPHI2FUN <- function(x, Beta) return(Beta[1] + Beta[2] * x)
	SIGMAFUN <- function(x, y) return((summary(lm(y ~ x))$sigma) ^ 2)
	
	# window #
	WindowLength <- length(x) # number of windows
	WindowCenter <- seq(1, WindowLength, by = step)
	WindowStart <- ((WindowCenter - width) <= 0) * 1 + 
		((WindowCenter - width) > 0) * 1 * (WindowCenter - width)
	WindowEnd <- ((WindowLength - WindowCenter) >= width) * 1 * (WindowCenter + width) +
		((WindowLength - WindowCenter) < width) * 1 * WindowLength
	x.list <- list()
	y.list <- list()
	for(i in seq(WindowLength)) {
		x.list[[i]] <- x[WindowStart[i] : WindowEnd[i]]
		y.list[[i]] <- y[WindowStart[i] : WindowEnd[i]]
	}
	
	# Bayesian #
	Linear.beta <- matrix(0, 2, WindowLength)
	for(i in seq(WindowLength)) {
		Linear.beta[ , i] <- BETAFUN(x.list[[i]], y.list[[i]])
	}
	Linear.sigma2 <- as.vector(mapply(SIGMAFUN, x = x.list, y = y.list))
#	for(i in seq(length(Linear.sigma2))){
#		if(is.nan(Linear.sigma2[i])) {
#			Linear.sigma2[i] <- runif(1, 0, 0.0001)
#		}
#	}
	Linear.beta.mean <- apply(Linear.beta, 1, mean)
#	for(i in seq(ncol(Linear.beta))){
#		if(sum(abs(Linear.beta.mean - Linear.beta[ , i])) == 0){
#			temp <- rnorm(length(Linear.beta.mean), 0, 0.001)
#			Linear.beta.mean <- Linear.beta.mean + temp
#		}
#	}
	A <- matrix(0, 2, 2)
	for(i in seq(WindowLength)) {
		A <- A + tcrossprod(Linear.beta[ , i] - Linear.beta.mean) / Linear.sigma2[i] 
	}
	V.beta <- A / WindowLength
	invgamma.mu <- mean(Linear.sigma2)
	invgamma.sigma <- sqrt(var(Linear.sigma2))
	a <- 2 + invgamma.mu / invgamma.sigma
	b <- (1 + invgamma.mu / invgamma.sigma) * invgamma.mu
	BayesNIG.beta <- matrix(0, 2, WindowLength)
	for(i in seq(WindowLength)) {
		x.t <- x.list[[i]]
		X.t <- cbind(rep(1, length(x.t)), x.t)
		y.t <- y.list[[i]]
		mu.star <- ginv(ginv(V.beta) + crossprod(X.t)) %*%
			(ginv(V.beta) %*% Linear.beta.mean + crossprod(X.t, y.t))
		BayesNIG.beta[ , i] <- mu.star
	}
	result <- list(Beta = BayesNIG.beta)
	return(result)
}

BayesianNIG <- function(x, ...) UseMethod("BayesianNIG")

BayesianNIG.Pheno <- function(object){
	
	data <- object$WorkingDataTemp
	step <- object$step
	width <- object$width
	label <- object$label
	x <- object$x 
	y <- object$y
	labelUniqueNames <- object$labelUniqueNames

	Beta <- NULL
	
	for(i in seq(length(x))){
		for(j in seq(length(y))){
			featuresTEMP <- c(x[i], y[j])
			if(!is.null(label)){
				splitDatabyLabels <- split(data, data[[label]])
				splitDataLabelNames <- names(splitDatabyLabels)
				interceptSplitData <- c()
				slopeSplitData <- c()
				labelSplitData <- c()
				for(k in seq(length(labelUniqueNames))){
					DataTemp <- splitDatabyLabels[[k]] 
					xTemp <- DataTemp[[featuresTEMP[1]]]
					yTemp <- DataTemp[[featuresTEMP[2]]]
					betaTemp <- BAYNIGFUN(xTemp, yTemp, step, width)$Beta
					interceptTemp <- as.vector(betaTemp[1, ])
					interceptSplitData <- c(interceptSplitData, interceptTemp)
					slopeTemp <- as.vector(betaTemp[2, ])
					slopeSplitData <- c(slopeSplitData, slopeTemp)
					labelSplitData <- c(labelSplitData, rep(splitDataLabelNames[k], 
						length(interceptTemp)))
				}
				BetaTemp <- cbind(interceptSplitData, slopeSplitData)
				colnames(BetaTemp) <- c(paste(x[i], y[j], "I", sep = "_"),
					paste(x[i], y[j], "S", sep = "_"))
				Beta <- cbind(Beta, BetaTemp)
			} else{
				DataTemp <- data 
				xTemp <- DataTemp[[featuresTEMP[1]]]
				yTemp <- DataTemp[[featuresTEMP[2]]]
				betaTemp <- BAYNIGFUN(xTemp, yTemp, step, width)$Beta
				interceptTemp <- as.vector(betaTemp[1, ])
				slopeTemp <- as.vector(betaTemp[2, ])
				BetaTemp <- cbind(interceptTemp, slopeTemp)
				colnames(BetaTemp) <- c(paste(x[i], y[j], "I", sep = "_"),
					paste(x[i], y[j], "S", sep = "_"))
				Beta <- cbind(Beta, BetaTemp)
			}
		}
	}
	Beta <- as.data.frame(Beta)	
	if(!is.null(label)){
		Beta <- data.frame(Beta, Label = labelSplitData)
	}
	return(Beta)
}	

BayesianNIG.predictPheno <- function(object){
	WorkingData <- object$data
	step <- object$step
	width <- object$width
#	label <- object$label
	x <- object$x 
	y <- object$y
#	labelUniqueNames <- object$labelUniqueNames
	
	Beta <- NULL
	for(i in seq(length(x))){
		for(j in seq(length(y))){
			featuresTEMP <- c(x[i], y[j])
			DataTemp <- WorkingData 
			xTemp <- DataTemp[[featuresTEMP[1]]]
			yTemp <- DataTemp[[featuresTEMP[2]]]
			betaTemp <- BAYNIGFUN(xTemp, yTemp, step, width)$Beta
			interceptTemp <- as.vector(betaTemp[1, ])
			slopeTemp <- as.vector(betaTemp[2, ])
			BetaTemp <- cbind(interceptTemp, slopeTemp)
			colnames(BetaTemp) <- c(paste(x[i], y[j], "I", sep = "_"),
				paste(x[i], y[j], "S", sep = "_"))
			Beta <- cbind(Beta, BetaTemp)
		}
	}
	Beta <- as.data.frame(Beta)

}

summary.Pheno <- function(object){
	
	topFeatureNames <- object$fea
	resPCA <- object$resPCA
	numberPCA <- object$numComp
	infTemp <- paste("Recommended individual feature(s) : ", 
		topFeatureNames, sep = " ")
	for(i in seq(length(infTemp))){
		print(infTemp[i])
	}
	infTemp <- paste("Number of selected component(s) in PCA : ", numberPCA, 
		sep = " ")
	print(infTemp)
	print("Information on PCA : ")
	return(summary(resPCA))
}

summary.cvPheno <- function(object){
	
	outputTable <- object$outputTable
	return(print(outputTable))
}

predict.Pheno <- function(object, data, x, y, 
	block, orderby, step = 1, width = 6){
		
	fn <- "predictPheno"
	arugmentsData <- list(data = data, x = x, y = y, block = block, 
		orderby = orderby, step = step, width = width)
	attr(arugmentsData, 'class') <- fn
	CheckArguments(arugmentsData)
	

#	orginalData <- object$org
#	BayesianData <- object$Bay
	predictData <- object$pre
#	clusterData <- object$clu
	topFeatureNames <- object$fea
	topFeatureFullNames <- object$feaf
	argumentsDataTemp <- object$arg
	labelTemp <- argumentsDataTemp$label
	
	if(is.null(labelTemp))
		stop("invalid usage, prediction is not available for clustering")
	methodTemp <- argumentsDataTemp$method
	xSelect <- intersect(x, topFeatureFullNames)
	ySelect <- intersect(y, topFeatureFullNames)
	
	valueTransferData <- TransferData(data, xSelect, ySelect, label = NULL, block, orderby)
	WorkingData <- valueTransferData$data
	
	objectlist <- list(data = WorkingData, x = xSelect, y = ySelect, step = step, width = width)
	attr(objectlist, 'class') <- fn
	WorkingDataTemp <- BayesianNIG(objectlist)
	
	inputData <- subset(WorkingDataTemp, select = topFeatureNames)
	
	if(methodTemp == "RF"){
		par.pred <- predict(predictData, inputData)	
	} else{
		par.pred <- predict(predictData, inputData)
	}
	outputData <- data.frame(inputData, Label = par.pred)
	
	return(outputData)
}

cv <- function(x, ...) UseMethod("cv")

cv.Pheno <- function(object, cvNumber, testBlockProp){
	
	# selecting feature data from original data
	# splitting
	
	data <- object$org
	WorkingData <- data
	argumentsDataTemp <- object$arg
	block <- argumentsDataTemp$block
	label <- argumentsDataTemp$label
	orderby <- argumentsDataTemp$orderby
	step <- argumentsDataTemp$step
	width <- argumentsDataTemp$width
	method <- argumentsDataTemp$method
	defaultLabel <- argumentsDataTemp$defaultLabel
	blockNamesONE <- argumentsDataTemp$blockNamesONE
	blockNamesTWO <- argumentsDataTemp$blockNamesTWO
	
	topFeatureFullNames <- object$feaf
	topFeatureNames <- object$fea
	x <- intersect(argumentsDataTemp$x, topFeatureFullNames)
	y <- intersect(argumentsDataTemp$y, topFeatureFullNames)
	
	blockUniqueNames <- unique(WorkingData[[block]])
	blockOrderedNames <- mixedsort(blockUniqueNames)
	labelUniqueNames <- unique(WorkingData[[label]])
	defaultLabelUniqueNames <- defaultLabel
	inverseLabel <- labelUniqueNames[-which(labelUniqueNames == defaultLabel)]
	
	splitDatabyBlock <- split(WorkingData, WorkingData[[block]])
	FUNALLEQU <- function(x, label){
		return(!any(x[[label]] != x[[label]][1]))
	}
	FUNLABEL <- function(x, label){
		return(x[[label]][1])
	}
	if(!all(sapply(splitDatabyBlock, FUNALLEQU, label = label, simplify = TRUE)))
		stop("data in some 'block' have multiple 'label'")
	blockOrderedLabels <- as.vector(sapply(splitDatabyBlock, FUNLABEL, 
		label = label, simplify = TRUE))
	
	WorkingDataSub <- subset(WorkingData, 
		select = c(x, y, label, block, orderby))
	features <- c(x, y)
	inicvNumber <- 1
	outputTable <- data.frame(matrix(0, length(blockUniqueNames), 3)) 
	rownames(outputTable) <- blockOrderedNames
	colnames(outputTable) <- c("performance", "precision", "recall")
	outputTableTemp <- outputTable 
	countTable <- rep(0, length(blockUniqueNames))
	while(inicvNumber <= cvNumber){
		cat("computing ", inicvNumber, "\r")
		valueSplitData <- SplitData(WorkingDataSub, block, orderby, 
			testBlockProp, blockOrderedLabels, blockOrderedNames)
		dataSplitData <- valueSplitData$data
		trainData <- dataSplitData$train
		testData <- dataSplitData$test
		testName <- dataSplitData$testName
		
		trainDataTemp <- subset(trainData, select = c(features, label))
		objectlist <- list(WorkingDataTemp = trainDataTemp, step = step,
			width = width, label = label, x = x, y = y,
			labelUniqueNames = labelUniqueNames)
		attr(objectlist, 'class') <- "Pheno"
		trainBayesianData <- BayesianNIG(objectlist)
		trainBayesianDataTemp <- subset(trainBayesianData, 
			select = c(topFeatureNames, "Label"))
		
		testDataTemp <- subset(testData, select = c(features, label))
		objectlist <- list(WorkingDataTemp = testDataTemp, step = step,
			width = width, label = label, x = x, y = y,
			labelUniqueNames = labelUniqueNames)
		attr(objectlist, 'class') <- "Pheno"
		testBayesianData <- BayesianNIG(objectlist)
		testBayesianDataTemp <- subset(testBayesianData, 
			select = c(topFeatureNames, "Label"))
		inputData <- subset(testBayesianDataTemp, select = topFeatureNames)
		
		if(method == "RF"){
			predictData <- randomForest(Label ~ ., data = trainBayesianDataTemp, 
				importance = TRUE, proximity = TRUE)
			par.pred <- predict(predictData, inputData)	
		} else{
			predictData <- svm(Label ~ ., data = trainBayesianDataTemp)
			par.pred <- predict(predictData, inputData)
		}
		
		DI <- DD <- II <- ID <- 0
		for(i in seq(length(par.pred))){
			if(par.pred[i] == testBayesianDataTemp$Label[i]){
				if(par.pred[i] %in% inverseLabel){
					II <- II + 1
				} else{
					DD <- DD + 1
				}
			} else{
				if(par.pred[i] == defaultLabel){
					DI <- DI + 1			
				} else{
					if(testBayesianDataTemp$Label[i] == defaultLabel){
						ID <- ID + 1
					}
				}
			}
		}	
		performance <- as.numeric((DD + II) / length(par.pred), 4)
		recall <- as.numeric(DD / length(which(testBayesianDataTemp$Label == defaultLabel)), 4) 
		precision <- as.numeric(DD / (DD + DI), 4)
		
		inicvNumber <- inicvNumber + 1
		rowID <- which(rownames(outputTable) %in% testName)
		countTable[rowID] <- countTable[rowID] + 1
		outputTableTemp[rowID, 1] <- outputTableTemp[rowID, 1] + 
			performance
		outputTableTemp[rowID, 2] <- outputTableTemp[rowID, 2] + 
			precision
		outputTableTemp[rowID, 3] <- outputTableTemp[rowID, 3] + 
			recall
	}
	
	outputTable <- outputTableTemp / countTable
	returnData <- list(
		outputTable = outputTable,
		blockNamesONE = blockNamesONE,
		blockNamesTWO = blockNamesTWO,
		topFeatureNames = topFeatureNames
		)
	attr(returnData,'class') <- "cvPheno"
	return(returnData)
}

plot.cvPheno <- function(object){
	
	blocknamesONE <- object$blockNamesONE
	blocknamesTWO <- object$blockNamesTWO
	if(is.null(blocknamesONE))
		stop("invalid 'block', which must have two components")
	topFeatureNames <- object$topFeatureNames
	output <- object$outputTable
	Lrow <- length(blocknamesONE)
	Lrun <- length(blocknamesTWO)
	
	rowN <- c()
	for(i in seq(Lrow)){
		rowN <- c(rowN, rep(blocknamesONE[i], Lrun))
	}
	colN <- rep(blocknamesTWO, Lrow)
	Values <- c(output[ , 1], output[ , 2], output[ , 3])
	Series <- c(rep("performance", Lrow * Lrun), rep("precision", Lrow * Lrun), 
		rep("recall", Lrow * Lrun))
	A <- data.frame(rowN, colN, Values, Series)
	avgPerformance <- round(mean(output[ , 1]), 4)
	avgPrecision <- round(mean(output[ , 2]), 4)
	avgRecall <- round(mean(output[ , 3]), 4)
	
	gg <- ggplot(A, aes(x = colN, y = rowN, fill = Values))
	gg <- gg + geom_tile(color = "white", size = 0.1)
	gg <- gg + scale_fill_viridis(name = "rate")
	gg <- gg + coord_equal()
	gg <- gg + facet_wrap( ~ Series, ncol = 1)
	gg <- gg + labs(x = NULL, y = NULL, 
		title = paste("Feature: ", toString(topFeatureNames), "\n", 
			"Averaged performance ", avgPerformance, 
			", precision ", avgPrecision, 
			", recall ", avgRecall))
	gg <- gg + theme_tufte(base_family = "sans")
	gg <- gg + theme(axis.ticks = element_blank())
	gg <- gg + theme(axis.text = element_text(size = 10))
	gg <- gg + theme(panel.border = element_blank())
	gg <- gg + theme(plot.title = element_text(hjust = 0))
	gg <- gg + theme(strip.text = element_text(size = 15, hjust = 0))
	gg <- gg + theme(panel.margin.x = unit(0.5, "cm"))
	gg <- gg + theme(panel.margin.y = unit(0.5, "cm"))
	gg <- gg + theme(legend.title = element_text(size = 10))
	gg <- gg + theme(legend.title.align = 1)
	gg <- gg + theme(legend.text = element_text(size = 10))
	gg <- gg + theme(legend.position = "bottom")
	gg <- gg + theme(legend.key.size = unit(0.2, "cm"))
	gg <- gg + theme(legend.key.width = unit(1, "cm"))
	gg
}

Pheno <- function(data = NULL, x = NULL, y = NULL, label = NULL, 
	defaultLabel = NULL, ncluster = NULL, block = NULL, orderby = NULL, 
	method = "SVM",	step = 1, width = 6, nfeatures = 3){
	
	fn <- "Pheno"
	arugmentsData <- list(data = data, x = x, y = y, label = label, 
		defaultLabel = defaultLabel, block = block, ncluster = ncluster,
		orderby = orderby, method = method, step = step, width = width,
		nfeatures = nfeatures)
	attr(arugmentsData, 'class') <- fn
	CheckArguments(arugmentsData)
	
	x <- unique(x)
	y <- unique(y)
	ret.x <- x
	ret.y <- y
	ret.label <- label
	ret.defaultLabel <- defaultLabel
	ret.block <- block
	ret.orderby <- orderby
	
	if(length(ret.block) == 2){
		blockNamesONE <- unique(data[[block[1]]])
		blockNamesTWO <- unique(data[[block[2]]])
	} else{
		blockNamesONE <- NULL
		blockNamesTWO <- NULL
	}
	
	valueTransferData <- TransferData(data, x, y, label, block, orderby)
	WorkingData <- valueTransferData$data
	orginalData <- WorkingData
	block <- valueTransferData$block
	orderby <- valueTransferData$orderby
	
	if(!is.null(label)){
		labelUniqueNames <- unique(WorkingData[[label]])
	
		features <- c(x, y)
		WorkingDataTemp <- subset(WorkingData, select = c(features, label))
#		orginalData <- WorkingDataTemp
#		crossvadalitionData <- WorkingData
		
		objectlist <- list(WorkingDataTemp = WorkingDataTemp, step = step,
			width = width, label = label, x = x, y = y,
			labelUniqueNames = labelUniqueNames)
		attr(objectlist, 'class') <- fn
		WorkingDataTemp <- BayesianNIG(objectlist)
		BayesianData <- WorkingDataTemp
		
		resPCA <- prcomp(subset(WorkingDataTemp, select = - Label), 
			center = TRUE, scale. = TRUE)
		propvar <- as.vector(summary(resPCA)$importance[2, ])
		cumuvar <- as.vector(summary(resPCA)$importance[3, ])
		numberPCA <- max(2, which(cumuvar > 0.7)[1])
		resRotation <- resPCA$rotation[ , (1 : (numberPCA))] 
		nameTemp <- rownames(resRotation)
		nameSplit <- strsplit(nameTemp, "_")
		scaleRotation <- t(t(abs(resRotation)) * propvar)
		
		# searching top feature names combining (I, S)
#		rowSumsTemp <- as.vector(rowSums(scaleRotation))
#		feaSums <- rep(0, length(rowSumsTemp) / 2)
#		feaNams <- rep(0, length(rowSumsTemp) / 2)
#		for(i in seq(length(rowSumsTemp) / 2)){
#			feaSums[i] <- rowSumsTemp[2 * i - 1] + rowSumsTemp[2 * i]
#			feaNams[i] <- paste(nameSplit[[2 * i]][1], nameSplit[[2 * i]][2], 
#				sep = "_")
#		}
#		featureTemp <- data.frame(feaNams, feaSums)
#		featureTemp <- featureTemp[order(featureTemp$feaSums, decreasing = TRUE), ]
#		topFeatureNames <- as.vector(featureTemp$feaNams[1 : min(nrow(featureTemp), 3)])
		
		# searching top feature names not combining (I, S)
		
		rowSumsTemp <- as.vector(rowSums(scaleRotation))
		feaSums <- rowSumsTemp
		feaNams <- as.vector(nameTemp)
		featureTemp <- data.frame(feaNams, feaSums)
		featureTemp <- featureTemp[order(featureTemp$feaSums, decreasing = TRUE), ]	
		topFeatureNames <- as.vector(featureTemp$feaNams[1 : min(nrow(featureTemp), nfeatures)])
		
		fullNameSplit <- strsplit(topFeatureNames, "_")
		topFeatureFullNames <- c()
		for(i in seq(length(topFeatureNames))){
			topFeatureFullNames <- c(fullNameSplit[[i]][1], fullNameSplit[[i]][2],
				topFeatureFullNames)
		}
		topFeatureFullNames <- unique(as.vector(topFeatureFullNames))
		
#		WorkingDataTempPCA <- as.matrix(WorkingDataTemp[ , - ncol(WorkingDataTemp)]) %*%
#			resRotation
#		WorkingDataTempPCA <- as.data.frame(WorkingDataTempPCA)	
#		WorkingDataTempPCA <- data.frame(WorkingDataTempPCA, Label = 
#			WorkingDataTemp[ , ncol(WorkingDataTemp)])
		WorkingDataTempPCA <- subset(BayesianData, select = c(topFeatureNames, 
			"Label"))
		
		if(method == "RF"){
			predictData <- randomForest(Label ~ ., data = WorkingDataTempPCA, 
				importance = TRUE, proximity = TRUE)
		} else{
			predictData <- svm(Label ~ ., data = WorkingDataTempPCA)
		}
		clusterData <- NULL
	} else{
		stop("invalid 'label' which is NULL")
#		features <- c(x, y)
#		labelUniqueNames <- NULL
#		WorkingDataTemp <- subset(WorkingData, select = c(features))
#		orginalData <- WorkingDataTemp
#		WorkingDataTemp <- BayesianNIGdirect(WorkingDataTemp, 
#			step, width, label, features, labelUniqueNames)
#		BayesianData <- WorkingDataTemp
#		predictData <- NULL
#		clusterResult <- kmeans(x = WorkingDataTemp, centers = ncluster)
#		membershipvector <- clusterResult$cluster
#		clusterData <- data.frame(BayesianData, Label = membershipvector)
	}
	
	returnData <- list(
		org = orginalData,
		Bay = BayesianData, # add topFeatureNames as their column names
		pre = predictData,
#		clu = clusterData,
#		cv = crossvadalitionData,
		arg = list(label = label, x = x, y = y, method = method,
			block = block, defaultLabel = defaultLabel, orderby = orderby,
			step = step, width = width, blockNamesONE = blockNamesONE,
		blockNamesTWO = blockNamesTWO),
#		pca = WorkingDataTempPCA, # matrix after rotation
		resPCA = resPCA,
		numComp = numberPCA,
		fea = topFeatureNames,
		feaf = topFeatureFullNames
		)
	attr(returnData,'class') <- fn
	return(returnData)
	
}

