CheckDataFrame <- function(data){
	ValueDataFrame <- is.data.frame(data)
	if(ValueDataFrame != TRUE){
		report <- "Error 01"
		returnData <- NULL
	} else{
		report <- "Successfully"
		returnData <- data
	}
	return(list(report = report, data = returnData))
}
CheckMissingValue <- function(data){ 
	if(nrow(na.omit(data)) == 0){
		report <- "Error 02"
		returnData <- NULL
	} else{
		report <- "Successfully"
		returnData <- na.omit(data)
	}
	return(list(report = report, data = returnData))
}
CheckFeature <- function(data, features, width){
	if(length(features) > 2){
		report <- "Error 04"
		returnData <- NULL
	} else{
		if(width > nrow(data)){
			report <- "Error 05: width is too large"
			returnData <- NULL
		} else{
			if(width < 4){
				report <- "Error 06: width is too small"
				returnData <- NULL
			} else{
				report <- "Successfully"
				returnData <- data
			}
		}
	}
	return(list(report = report, data = returnData))
}
CheckBlock <- function(data, block, label){
	if(is.null(block) == TRUE){
		report <- "Error 07"
		returnData <- NULL
		returnLabelName <- NULL
	} else{
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
			report <- "Successfully"
			returnData <- data
			returnLabelName <- blockLabelName
		} else{
			report <- "Error block has multiple labels"
			returnData <- NULL
			returnLabelName <- NULL
		}
	}
	return(list(report = report, data = returnData, blockLabelName = returnLabelName))
}
SplitData <- function(data, blockLabelName, block, label, orderby, testBlockProp){
	if(testBlockProp > 1 | testBlockProp < 0){
		report <- "Error 08"
		returnData <- NULL
	} else{
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
		report <- "Successfully"
		returnData <- list(test = testData, train = trainData, 
			testName = testName)	
	}
	
	return(list(report = report, data = returnData))
}
BAYNIGFUN <- function(x, y, width) {
	# compute Bayesian NIG for any input without spliting
	require(MASS)
	R2FUN <- function(y, Est.y) return(1 - sum((y - Est.y) ^ 2) / sum((y - mean(y)) ^ 2))
	BETAFUN <- function(x, y) return(lm(y ~ x)$coefficients)
	BETAPHI2FUN <- function(x, Beta) return(Beta[1] + Beta[2] * x)
	SIGMAFUN <- function(x, y) return((summary(lm(y ~ x))$sigma) ^ 2)
	
	# window #
	WindowLength <- length(x) # number of windows
	WindowCenter <- seq(WindowLength)
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
	for(i in seq(length(Linear.sigma2))){
		if(is.nan(Linear.sigma2[i])) {
			Linear.sigma2[i] <- runif(1, 0, 0.0001)
		}
	}
	Linear.beta.mean <- apply(Linear.beta, 1, mean)
	for(i in seq(ncol(Linear.beta))){
		if(sum(abs(Linear.beta.mean - Linear.beta[ , i])) == 0){
			temp <- rnorm(length(Linear.beta.mean), 0, 0.001)
			Linear.beta.mean <- Linear.beta.mean + temp
		}
	}
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
BayesianNIG <- function(data, width, features, label) {
	trainData <- data$train
	testData <- data$test
	
	trainDataONE <- trainData[which(trainData[[label]] == 
		unique(trainData[[label]])[1]), ]
	xTrainONE <- trainDataONE[[features[1]]]
	yTrainONE <- trainDataONE[[features[2]]]
	BetaONE <- BAYNIGFUN(xTrainONE, yTrainONE, width)$Beta
	interceptTrainONE <- as.vector(BetaONE[1, ])
	slopeTrainONE <- as.vector(BetaONE[2, ])
	labelTrainONE <- rep(unique(trainData[[label]])[1], length(interceptTrainONE))
	
	trainDataTWO <- trainData[which(trainData[[label]] == 
		unique(trainData[[label]])[2]), ]
	xTrainTWO <- trainDataTWO[[features[1]]]
	yTrainTWO <- trainDataTWO[[features[2]]]
	BetaTWO <- BAYNIGFUN(xTrainTWO, yTrainTWO, width)$Beta
	interceptTrainTWO <- as.vector(BetaTWO[1, ])
	slopeTrainTWO <- as.vector(BetaTWO[2, ])
	labelTrainTWO <- rep(unique(trainData[[label]])[2], length(interceptTrainTWO))
	
	BetaTrain <- list(intercept = c(interceptTrainONE, interceptTrainTWO), 
		slope = c(slopeTrainONE, slopeTrainTWO), label = c(labelTrainONE, labelTrainTWO))
	
	testDataONE <- testData[which(testData[[label]] == 
		unique(testData[[label]])[1]), ]	
	xTestONE <- testDataONE[[features[1]]]
	yTestONE <- testDataONE[[features[2]]]
	BetaTest <- BAYNIGFUN(xTestONE, yTestONE, width)$Beta
	interceptTestONE <- as.vector(BetaTest[1, ])
	slopeTestONE <- as.vector(BetaTest[2, ])
	labelTestONE <- rep(unique(testData[[label]])[1], length(interceptTestONE))
	
	
	testDataTWO <- testData[which(testData[[label]] == 
		unique(testData[[label]])[2]), ]	
	xTestTWO <- testDataTWO[[features[1]]]
	yTestTWO <- testDataTWO[[features[2]]]
	BetaTest <- BAYNIGFUN(xTestTWO, yTestTWO, width)$Beta
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
	if(is.null(defaultLabel) == TRUE){
		report <- "Error 09"
		returnData <- NULL
	} else{
		if(defaultLabel == unique(InputTrain$label)[1]){
			defaultLabel <- unique(InputTrain$label)[1]
			inverseLabel <- unique(InputTrain$label)[2]
			report <- "Successfully"
		} else{
			if(defaultLabel == unique(InputTrain$label)[2]){
				defaultLabel <- unique(InputTrain$label)[2]
				inverseLabel <- unique(InputTrain$label)[1]
				report <- "Successfully"
			} else{
				report <- "Error 10"
				returnData <- NULL
			}
		}
		if(report != "Successfully"){
			report <- report
			returnData <- returnData
		} else{
			InputTrain <- data.frame(
				intercept = InputTrain$intercept,
				slope = InputTrain$slope, 
				label = InputTrain$label)
			InputTest <- data.frame(
				intercept = InputTest$intercept,
				slope = InputTest$slope, 
				label = InputTest$label)		
			if(method == "RF"){
				require(randomForest)
				par.rf <- randomForest(label ~ .,
					data = InputTrain, 
					importance = TRUE, 
					proximity = TRUE)
				par.pred <- predict(par.rf, subset(InputTest, select = - label))
			
			} else{
				if(method == "SVM"){
					require(e1071)
					par.svm <- svm(label ~ ., 
					data = InputTrain)
					par.pred <- predict(par.svm, subset(InputTest, select = - label))
				}
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
			report <- "Successfully"
			returnData <- list(
				precision = round(as.numeric(precision), 4),
				recall = round(as.numeric(recall), 4),
				performance = round(as.numeric(performance), 4))
		}
	}
	return(list(report = report, data = returnData))
}
PhenoPro <- function(data, label, defaultLabel,
	features, block, method,
	orderby, width, cvNumber, testBlockProp, visualization){
	
	valueCheckDataFrame <- CheckDataFrame(data) # check data frame
	if(valueCheckDataFrame$report != "Successfully"){
		return(valueCheckDataFrame$report)
	} else {
		dataCheckDataFrame <- valueCheckDataFrame$data
	}
	WorkingData <- subset(dataCheckDataFrame, 
		select = c(label, features, block, orderby))
		
	valueCheckMissingValue <- CheckMissingValue(WorkingData) # check NA
	if(valueCheckMissingValue$report != "Successfully"){
		return(valueCheckMissingValue$report)
	} else{
		dataCheckMissingValue <- valueCheckMissingValue$data
	}
	WorkingData <- dataCheckMissingValue
	
	WorkingData <- WorkingData[order(WorkingData[[orderby]], decreasing = FALSE), ]
	
	valueCheckFeature <- CheckFeature(WorkingData, features, width) # check features
	if(valueCheckFeature$report != "Successfully"){
		return(valueCheckFeature$report)
	} else{
		dataCheckFeature <- valueCheckFeature$data
	}
	WorkingData <- dataCheckFeature
	
	valueCheckBlock <- CheckBlock(WorkingData, block, label)
	if(valueCheckBlock$report != "Successfully"){
		return(valueCheckBlock$report)
	} else{
		dataCheckBlock <- valueCheckBlock$data
		blockLabelName <- valueCheckBlock$blockLabelName
	}
	WorkingData <- dataCheckBlock
	
	outputTable <- data.frame(matrix(0, length(unique(WorkingData[[block]])), 3))
	countTabel <- rep(0, length(unique(WorkingData[[block]])))
	rownames(outputTable) <- unique(WorkingData[[block]])[order(unique(WorkingData[[block]]))]
	colnames(outputTable) <- c("performance", "precision", "recall")
	inicvNumber <- 1
	while(inicvNumber <= cvNumber){
		cat('cross validation number:', inicvNumber, "\r")
		valueSplitData <- SplitData(WorkingData, blockLabelName, block, 
			label, orderby, testBlockProp)
		if(valueSplitData$report != "Successfully"){
			inicvNumber <- inicvNumber
		} else{
			dataSplitData <- valueSplitData$data
			WorkingDataTemp <- dataSplitData
			testNameTemp <- WorkingDataTemp$testName
		}
		
		if(is.null(label) != TRUE){
			WorkingDataTemp <- BayesianNIG(WorkingDataTemp, width, features, label)
			valuePrediction <- Prediction(WorkingDataTemp, method, label, defaultLabel)
			if(valuePrediction$report != "Successfully"){
				inicvNumber <- inicvNumber
			} else{
				dataPrediction <- valuePrediction$data
				testName <- testNameTemp
				inicvNumber <- inicvNumber + 1
				WorkingDataTemp <- dataPrediction
				rowID <- which(rownames(outputTable) %in% testName)
				countTabel[rowID] <- countTabel[rowID] + 1
				outputTable[rowID, 1] <- outputTable[rowID, 1] + 
					WorkingDataTemp$performance
				outputTable[rowID, 2] <- outputTable[rowID, 2] + 
					WorkingDataTemp$precision
				outputTable[rowID, 3] <- outputTable[rowID, 3] + 
					WorkingDataTemp$recall
			}	
		}
	}
	output <- outputTable / countTabel
	return(list(output = output, outputTable = outputTable, countTabel = countTabel))
}



