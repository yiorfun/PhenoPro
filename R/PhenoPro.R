CheckInputArguments <- function(data, x, y, label, defaultLabel, block, orderby, 
	method, step, width, cvNumber, testBlockProp, visualization){
	
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
	if(!is.null(label) & length(unique(data[[label]])) < 2)
		stop("invalid 'label', which at least has two different components")
	if(!is.null(defaultLabel) & length(defaultLabel) >= length(unique(data[[label]])))
		stop("invalid 'defaultLabel', it must have less components than 'label'")	
	if(!is.null(defaultLabel) & !(defaultLabel %in% unique(data[[label]])))
		stop("invalid 'defaultLabel', which should be any component of 'label'")
	
	# check block argument
	if(!is.null(block) & length(block) > 2)
		stop("invalid 'block', at most 2 dimensions")
	
	# check orderby argument
	if(!is.null(orderby) & length(orderby) >= 2)
		stop("invalid 'orderby', at most 1 dimension")
	
	# check method argument	
	if(!(method %in% c("SVM", "RF")))
		stop("invalid 'method', only 'SVM' and 'RF' are available")
	if(length(method) != 1)
		stop("invalid 'method', at most 1 dimension")
	
	# check step and width arguments
	if(!(step == round(step)))
		stop("invalid 'step', it must be an integer")
	if(!(width == round(width)))
		stop("invalid 'width', it must be an integer")
	
	# check cvNumber and testBlockProp arguments
	if(!(cvNumber == round(cvNumber)))
		stop("invalid 'cvNumber', it must be an integer")
	if(testBlockProp < 0 | testBlockProp > 1)
		stop("invalid 'testBlockProp', it must be between 0 and 1")
	
	# check visualization argument
	if(!is.logical(visualization))
		stop("invalid 'visualization', it must be a logical value (TorF)")
}

TransferData <- function(data, x, y, label, defaultLabel, block, orderby){
	
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
	if(!is.null(label) & !is.null(defaultLabel)){
		WorkingDataSub <- subset(WorkingData, select = c(
			x, y, label, block, orderby))
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
	testSample <- mapply(FUNSAMPLE, x = countLabelsTemp, y = testSizes)
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

BayesianNIGcv <- function(data, step, width, features, label,
	labelUniqueNames) {
	trainData <- data$train
	testData <- data$test

	splitTrainDatabyLabels <- split(trainData, trainData[[label]])
	splitTrainDataLabelNames <- names(splitTrainDatabyLabels)
	splitTestDatabyLabels <- split(testData, testData[[label]])
	splitTestDataLabelNames <- names(splitTestDatabyLabels)
	interceptTrain <- c()
	slopeTrain <- c()
	labelTrain <- c()
	interceptTest <- c()
	slopeTest <- c()
	labelTest <- c()
	for(i in seq(length(labelUniqueNames))){
		trainDataTemp <- splitTrainDatabyLabels[[i]] 
		xTrainTemp <- trainDataTemp[[features[1]]]
		yTrainTemp <- trainDataTemp[[features[2]]]
		betaTrainTemp <- BAYNIGFUN(xTrainTemp, yTrainTemp, step, width)$Beta
		interceptTrainTemp <- as.vector(betaTrainTemp[1, ])
		interceptTrain <- c(interceptTrain, interceptTrainTemp)
		slopeTrainTemp <- as.vector(betaTrainTemp[2, ])
		slopeTrain <- c(slopeTrain, slopeTrainTemp)
		labelTrain <- c(labelTrain, rep(splitTrainDataLabelNames[i], 
			length(interceptTrainTemp)))
			
		testDataTemp <- splitTestDatabyLabels[[i]] 
		xTestTemp <- testDataTemp[[features[1]]]
		yTestTemp <- testDataTemp[[features[2]]]
		betaTestTemp <- BAYNIGFUN(xTestTemp, yTestTemp, step, width)$Beta
		interceptTestTemp <- as.vector(betaTestTemp[1, ])
		interceptTest <- c(interceptTest, interceptTestTemp)
		slopeTestTemp <- as.vector(betaTestTemp[2, ])
		slopeTest <- c(slopeTest, slopeTestTemp)
		labelTest <- c(labelTest, rep(splitTestDataLabelNames[i], 
			length(interceptTestTemp)))
	}
	betaTrain <- list(intercept = interceptTrain, 
		slope = slopeTrain, 
		label = labelTrain)
	betaTest <- list(intercept = interceptTest, 
		slope = slopeTest, 
		label = labelTest)
	
	returnData <- list(train = betaTrain, test = betaTest)
	return(returnData)
}

Prediction <- function(data, method, label, defaultLabel, inverseLabel){
	InputTrain <- data$train
	InputTest <- data$test

	InputTrain <- data.frame(
		intercept = InputTrain$intercept,
		slope = InputTrain$slope, 
		label = InputTrain$label)
	InputTest <- data.frame(
		intercept = InputTest$intercept,
		slope = InputTest$slope, 
		label = InputTest$label)		
	if(method == "RF"){
#		require(randomForest)
		par.rf <- randomForest(label ~ ., data = InputTrain, 
			importance = TRUE, proximity = TRUE)
		par.pred <- predict(par.rf, subset(InputTest, select = - label))	
		} else{
#			require(e1071)
			par.svm <- svm(label ~ ., data = InputTrain)
			par.pred <- predict(par.svm, subset(InputTest, select = - label))
		}
	DI <- DD <- II <- ID <- 0
	for(i in seq(length(par.pred))){
		if(par.pred[i] == InputTest$label[i]){
			if(par.pred[i] %in% inverseLabel){
				II <- II + 1
			} else{
				DD <- DD + 1
			}
		} else{
			if(par.pred[i] == defaultLabel){
				DI <- DI + 1			
			} else{
				if(InputTest$label[i] == defaultLabel){
					ID <- ID + 1
				}
			}
		}
	}
	performance <- (DD + II) / length(par.pred)
	recall <- DD / length(which(InputTest$label == defaultLabel)) 
	precision <- DD / (DD + DI)
	returnData <- list(
		precision = round(as.numeric(precision), 4),
		recall = round(as.numeric(recall), 4),
		performance = round(as.numeric(performance), 4))
	return(list(data = returnData))
		
}


BayesianNIGdirect <- function(data, step, width, label, features, 
	labelUniqueNames){
	
	splitDatabyLabels <- split(data, data[[label]])
	splitDataLabelNames <- names(splitDatabyLabels)
	interceptSplitData <- c()
	slopeSplitData <- c()
	labelSplitData <- c()

	for(i in seq(length(labelUniqueNames))){
		DataTemp <- splitDatabyLabels[[i]] 
		xTemp <- DataTemp[[features[1]]]
		yTemp <- DataTemp[[features[2]]]
		betaTemp <- BAYNIGFUN(xTemp, yTemp, step, width)$Beta
		interceptTemp <- as.vector(betaTemp[1, ])
		interceptSplitData <- c(interceptSplitData, interceptTemp)
		slopeTemp <- as.vector(betaTemp[2, ])
		slopeSplitData <- c(slopeSplitData, slopeTemp)
		labelSplitData <- c(labelSplitData, rep(splitDataLabelNames[i], 
			length(interceptTemp)))
	}
	Beta <- data.frame(intercept = interceptSplitData, 
		slope = slopeSplitData, 
		label = labelSplitData)
	return(Beta)
}	


PhenoPro <- function(
	data = NULL, 
	x = NULL, 
	y = NULL, 
	label = NULL, 
	defaultLabel = NULL, 
	block = NULL, 
	orderby = NULL,
	method = "SVM",	
	step = 1, width = 5, 
	cvNumber = 100, 
	testBlockProp = 0.2, 
	visualization = FALSE){
	
	CheckInputArguments(data, x, y, label, defaultLabel, block, orderby,
		method,	step, width, cvNumber, testBlockProp, visualization)
	
	ret.x <- x
	ret.y <- y
	ret.label <- label
	ret.defaultLabel <- defaultLabel
	ret.block <- block
	ret.orderby <- orderby
	
	valueTransferData <- TransferData(data, x, y, label, defaultLabel, 
		block, orderby)
	WorkingData <- valueTransferData$data
	block <- valueTransferData$block
	orderby <- valueTransferData$orderby
	
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
	
	xcol <- length(x)
	ycol <- length(y)
	outputTable <- data.frame(matrix(0, length(blockUniqueNames), 3)) 
	
	rownames(outputTable) <- blockOrderedNames
	colnames(outputTable) <- c("performance", "precision", "recall")	
	output <- list()
	xycol <- 1
	for(i in seq(xcol)){
		for(j in seq(ycol)){
			outputTableTemp <- outputTable 
			countTable <- rep(0, length(blockUniqueNames))
			inicvNumber <- 1
			WorkingDataSub <- subset(WorkingData, 
				select = c(x[i], y[j], label, block, orderby))
			features <- c(x[i], y[j])
			outputNamesTemp <- paste(x[i], y[j], sep = "_") 
			while(inicvNumber <= cvNumber){
				cat(paste("computing ", outputNamesTemp), inicvNumber, "\r")
				valueSplitData <- SplitData(WorkingDataSub, block, orderby, 
	testBlockProp, blockOrderedLabels, blockOrderedNames)
				dataSplitData <- valueSplitData$data
				WorkingDataTemp <- dataSplitData
				testNameTemp <- WorkingDataTemp$testName
				
				WorkingDataTemp <- BayesianNIGcv(WorkingDataTemp, step, 
					width, features, label, labelUniqueNames)
				valuePrediction <- Prediction(WorkingDataTemp, method, 
					label, defaultLabel, inverseLabel)
				
				dataPrediction <- valuePrediction$data
				testName <- testNameTemp
				inicvNumber <- inicvNumber + 1
				WorkingDataTemp <- dataPrediction
				rowID <- which(rownames(outputTable) %in% testName)
				countTable[rowID] <- countTable[rowID] + 1
				outputTableTemp[rowID, 1] <- outputTableTemp[rowID, 1] + 
					WorkingDataTemp$performance
				outputTableTemp[rowID, 2] <- outputTableTemp[rowID, 2] + 
					WorkingDataTemp$precision
				outputTableTemp[rowID, 3] <- outputTableTemp[rowID, 3] + 
					WorkingDataTemp$recall
			}
			output[[outputNamesTemp]] <- outputTableTemp / countTable
			if(visualization){
				WorkingDataPlot <- subset(WorkingData, 
					select = c(x[i], y[j], label))
				PredictDataPlot <- BayesianNIGdirect(WorkingDataPlot, 
					step, width, label, features,
					labelUniqueNames)
				p1 <- ggplot() + 
					geom_point(data = WorkingDataPlot, 
					mapping = aes(WorkingData[[x[j]]], WorkingData[[y[i]]], 
						colour = factor(WorkingData[[label]]))) +
					xlab(x[i]) +
					ylab(y[j]) +
					theme(legend.position = "none")
				p2 <- ggplot() + 
					geom_point(data = PredictDataPlot, 
					mapping = aes(intercept, slope, 
						colour = label))
				p12 <- grid.arrange(p1, p2, ncol = 2)
				p12
			}
			xycol <- xycol + 1
		}
	}
	return(list(output = output))
}