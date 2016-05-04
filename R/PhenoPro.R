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
	if(!is.null(defaultLabel) & length(defaultLabel) >= 2)
		stop("invalid 'defaultLabel', at most 1 dimension")
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

CheckBlockStatus <- function(data, block, label){
	n <- length(unique(data[[block]]))
		blockOrderedNames <- unique(data[[block]])[order(unique(data[[block]]))]
		blockLabelNumber <- rep(0, n)
		blockLabelName <- rep(0, n)
		for(i in seq(n)){
			subdata <- data[which(data[[block]] == blockOrderedNames[i]) , ]
			blockLabelNumber[i] <- length(unique(subdata[[label]]))
			if(blockLabelNumber[i] == 1){
				blockLabelName[i] <- unique(subdata[[label]]) 
			} else{
				blockLabelName[i] <- NULL
			}
		}
		if(sum(blockLabelNumber) == n){
			returnData <- data
			returnLabelName <- blockLabelName
		} else{
			stop("'block' has multiple labels")
		}
	return(list(data = returnData, blockLabelName = returnLabelName))
}

SplitData <- function(data, blockLabelName, block, label, orderby, testBlockProp){
	blockOrderedNames <- unique(data[[block]])[order(unique(data[[block]]))]
	blockID <- seq(length(blockOrderedNames))
	blockDataFrame <- data.frame(Name = blockOrderedNames, ID = blockID, 
		Label = blockLabelName)
	subONE <- subset(blockDataFrame, Label == unique(blockLabelName)[1])
	subTWO <- subset(blockDataFrame, Label == unique(blockLabelName)[2])
	blockNumberONE <- nrow(subONE)
	blockNumberTWO <- nrow(subTWO)
	testSizeONE <- round(testBlockProp * blockNumberONE)
	testSizeTWO <- round(testBlockProp * blockNumberTWO)	
	testIDONE <- sample(blockNumberONE, size = testSizeONE, 
		replace = FALSE, prob = rep(1 / blockNumberONE, blockNumberONE))
	testIDONE <- testIDONE[order(testIDONE)]
	trainIDONE <- seq(blockNumberONE)[- testIDONE]
	testIDTWO <- sample(blockNumberTWO, size = testSizeTWO, 
		replace = FALSE, prob = rep(1 / blockNumberTWO, blockNumberTWO))
	testIDTWO <- testIDTWO[order(testIDTWO)]
	trainIDTWO <- seq(blockNumberTWO)[- testIDTWO]
		
	testName <- c(subONE[testIDONE, ]$Name, subTWO[testIDTWO, ]$Name)
	testName <- testName[order(testName)]
	testName <- blockOrderedNames[testName]
	trainName <- c(subONE[trainIDONE, ]$Name, subTWO[trainIDTWO, ]$Name)
	trainName <- trainName[order(trainName)]
	trainName <- blockOrderedNames[trainName]
		
	testData <- data[which(data[[block]] %in% testName), ]
	testData <- testData[order(testData[[orderby]]), ]
	trainData <- data[which(data[[block]] %in% trainName), ]
	trainData <- trainData[order(trainData[[orderby]]), ]
	returnData <- list(test = testData, train = trainData, 
		testName = testName)	
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

BayesianNIGcv <- function(data, step, width, features, label) {
	trainData <- data$train
	testData <- data$test
	
	trainDataONE <- trainData[which(trainData[[label]] == 
		unique(trainData[[label]])[1]), ]
	xTrainONE <- trainDataONE[[features[1]]]
	yTrainONE <- trainDataONE[[features[2]]]
	BetaONE <- BAYNIGFUN(xTrainONE, yTrainONE, step, width)$Beta
	interceptTrainONE <- as.vector(BetaONE[1, ])
	slopeTrainONE <- as.vector(BetaONE[2, ])
	labelTrainONE <- rep(unique(trainData[[label]])[1], length(interceptTrainONE))
	
	trainDataTWO <- trainData[which(trainData[[label]] == 
		unique(trainData[[label]])[2]), ]
	xTrainTWO <- trainDataTWO[[features[1]]]
	yTrainTWO <- trainDataTWO[[features[2]]]
	BetaTWO <- BAYNIGFUN(xTrainTWO, yTrainTWO, step, width)$Beta
	interceptTrainTWO <- as.vector(BetaTWO[1, ])
	slopeTrainTWO <- as.vector(BetaTWO[2, ])
	labelTrainTWO <- rep(unique(trainData[[label]])[2], length(interceptTrainTWO))
	
	BetaTrain <- list(intercept = c(interceptTrainONE, interceptTrainTWO), 
		slope = c(slopeTrainONE, slopeTrainTWO), label = c(labelTrainONE, labelTrainTWO))
	
	testDataONE <- testData[which(testData[[label]] == 
		unique(testData[[label]])[1]), ]	
	xTestONE <- testDataONE[[features[1]]]
	yTestONE <- testDataONE[[features[2]]]
	BetaTest <- BAYNIGFUN(xTestONE, yTestONE, step, width)$Beta
	interceptTestONE <- as.vector(BetaTest[1, ])
	slopeTestONE <- as.vector(BetaTest[2, ])
	labelTestONE <- rep(unique(testData[[label]])[1], length(interceptTestONE))
	
	testDataTWO <- testData[which(testData[[label]] == 
		unique(testData[[label]])[2]), ]	
	xTestTWO <- testDataTWO[[features[1]]]
	yTestTWO <- testDataTWO[[features[2]]]
	BetaTest <- BAYNIGFUN(xTestTWO, yTestTWO, step, width)$Beta
	interceptTestTWO <- as.vector(BetaTest[1, ])
	slopeTestTWO <- as.vector(BetaTest[2, ])
	labelTestTWO <- rep(unique(testData[[label]])[2], length(interceptTestTWO))
	
	BetaTest <- list(intercept = c(interceptTestONE, interceptTestTWO), 
		slope = c(slopeTestONE, slopeTestTWO), label = c(labelTestONE, labelTestTWO))
		
	returnData <- list(train = BetaTrain, test = BetaTest)
	return(returnData)
}

Prediction <- function(data, method, label, defaultLabel){
	InputTrain <- data$train
	InputTest <- data$test
	if(defaultLabel == unique(InputTrain$label)[1]){
		inverseLabel <- unique(InputTrain$label)[2]
	} else{
		inverseLabel <- unique(InputTrain$label)[1]
	}
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
			if(par.pred[i] == inverseLabel){
				II <- II + 1
			} else{
				DD <- DD + 1
			}
		} else{
			if(par.pred[i] == inverseLabel){
				ID <- ID + 1
			} else{
				DI <- DI + 1
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


BayesianNIGdirect <- function(data, x, y, step, width, label){
	dataONE <- data[which(data[[label]] == unique(data[[label]])[1]), ]
	xONE <- dataONE[[x]]
	yONE <- dataONE[[y]]
	BetaONE <- BAYNIGFUN(xONE, yONE, step, width)$Beta	
	interceptONE <- as.vector(BetaONE[1, ])
	slopeONE <- as.vector(BetaONE[2, ])
	labelONE <- rep(unique(data[[label]])[1], length(interceptONE))
	
	dataTWO <- data[which(data[[label]] == unique(data[[label]])[2]), ]
	xTWO <- dataTWO[[x]]
	yTWO <- dataTWO[[y]]
	BetaTWO <- BAYNIGFUN(xTWO, yTWO, step, width)$Beta
	interceptTWO <- as.vector(BetaTWO[1, ])
	slopeTWO <- as.vector(BetaTWO[2, ])
	labelTWO <- rep(unique(data[[label]])[2], length(interceptTWO))
	
	Beta <- data.frame(intercept = c(interceptONE, interceptTWO), 
		slope = c(slopeONE, slopeTWO), label = c(labelONE, labelTWO))
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
	
	x.raw <- x
	y.raw <- y
	label.raw <- label
	defaultLabel.raw <- defaultLabel
	block.raw <- block
	orderby.raw <- orderby
	
	valueTransferData <- TransferData(data, x, y, label, defaultLabel, block, orderby)
	WorkingData <- valueTransferData$data
	block <- valueTransferData$block
	orderby <- valueTransferData$orderby
	
	valueCheckBlockStatus <- CheckBlockStatus(WorkingData, block, label)
	blockLabelName <- valueCheckBlockStatus$blockLabelName
	WorkingData <- valueCheckBlockStatus$data
	
	xcol <- length(x)
	ycol <- length(y)
	outputTable <- data.frame(matrix(0, length(unique(WorkingData[[block]])), 3)) 
	
	rownames(outputTable) <- mixedsort(unique(WorkingData[[block]]))
	colnames(outputTable) <- c("performance", "precision", "recall")	
	output <- list()
	xycol <- 1
	for(i in seq(xcol)){
		for(j in seq(ycol)){
			outputTableTemp <- outputTable 
			countTable <- rep(0, length(unique(WorkingData[[block]])))
			inicvNumber <- 1
			WorkingDataSub <- subset(WorkingData, 
				select = c(x[i], y[j], label, block, orderby))
			features <- c(x[i], y[j])
			outputNamesTemp <- paste(x[i], y[j], sep = "_") 
			while(inicvNumber <= cvNumber){
				cat(paste("computing ", outputNamesTemp), inicvNumber, "\r")
				valueSplitData <- SplitData(WorkingDataSub, blockLabelName, block, 
					label, orderby, testBlockProp)
				dataSplitData <- valueSplitData$data
				WorkingDataTemp <- dataSplitData
				testNameTemp <- WorkingDataTemp$testName
				
				WorkingDataTemp <- BayesianNIGcv(WorkingDataTemp, step, width, features, label)
				valuePrediction <- Prediction(WorkingDataTemp, method, label, defaultLabel)
				
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
					x = x[i], y = y[j], step, width, label)
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

