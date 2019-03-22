
#' Function that obtains, from the genotype data, a dataset with the data
#' prepared for the Machine Learning part, after doing the variable selection
#'
#' @param workPath String with your work path
#' @param path2Geno Path to your genotype data
#' @param path2Covs Path to your covariates data
#' @param predictor What you want to predict
#' @param path2GWAS Path to your GWAS
#' @param path2PRSICE Path to PRSICE
#'
#' @return A handler with the mldata
#' @export
#'
#' @examples
#' handlerMLdata <- fromGenoToMLdata("/home/rafael",
#'                                   "/home/users/gsit/juanbot/JUAN_SpanishGWAS/UNRELATED.SPAIN4.HARDCALLS.Rsq0.8",
#'                                   "/home/users/gsit/juanbot/JUAN_SpanishGWAS/COVS_SPAIN",
#'                                   "DISEASE",
#'                                   "/home/users/gsit/juanbot/JUAN_SpanishGWAS/toJuanNov7th2018/",
#'                                   "/home/users/gsit/juanbot/genoml-core/otherPackages/")
fromGenoToMLdata = function(workPath,
                            path2Geno = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/UNRELATED.SPAIN4.HARDCALLS.Rsq0.8",
                            path2Covs = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/COVS_SPAIN",
                            predictor = "DISEASE",
                            path2GWAS = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/toJuanNov7th2018/",
                            path2PRSICE = "/home/users/gsit/juanbot/genoml-core/otherPackages/"
                            ){

  #Our first handler will always start by holding the genotype file, the covariates, the id and familial id columns
  #at the covariates, the variable name from which to predict how we want the phenotype file to be started
  h = getHandlerToGenotypeData(geno=path2Geno,
                               covs=path2Covs,
                               id="IID",
                               fid="FID",
                               predictor=predictor,
                               #With this we assure everything will be written under workPath
                               pheno=paste0(workPath,"/MyPhenotype"))

  # #Generate a holdout partition. The training part will be used to generate three ML models
  # #The test part will be used to evaluate the models (the handler is saved for later, not used in this script)
  holdout = getPartitionsFromHandler(genoHandler=h,
                                     workPath = workPath,
                                     path2plink="",
                                     how="holdout",
                                     p=0.50)

  holdout = genDataFromHandler(holdout,lazy=T)
  saveRDS(holdout,paste0(workPath,"/holdout.rds"))
  holdout = readRDS(paste0(workPath,"/holdout.rds"))

  handlerSNPs = mostRelevantSNPs(handler=getHandlerFromFold(handler=holdout,type="train",index=1),
                                 path2plink="",
                                 gwas="RISK_noSpain.tab",
                                 gwasDef=" --beta --snp MarkerName --A1 Allele1 --A2 Allele2 --stat Effect --se StdErr --pvalue P-value",
                                 path2GWAS=path2GWAS,
                                 PRSiceexe="PRSice_linux",
                                 path2PRSice=path2PRSice,
                                 clumpField = "P-value",
                                 SNPcolumnatGWAS = "MarkerName")

  saveRDS(handlerSNPs,paste0(workPath,"/handlerSNPs.rds"))
  handlerSNPs = readRDS(paste0(workPath,"/handlerSNPs.rds"))

  mldatahandler = fromSNPs2MLdata(handler=holdout,
                                  addit="NA",
                                  path2plink="",
                                  fsHandler=handlerSNPs)

  saveRDS(mldatahandler,paste0(workPath,"/mldatahandler.rds"))
  mldatahandler = readRDS(paste0(workPath,"/mldatahandler.rds"))

  #These are the two ML datasets generated from the genotype data and the
  #feature selection
  cat("Your train dataset is at",mldatahandler$train1mldata,"\n")
  cat("Your test dataset is at",mldatahandler$test1mldata,"\n")

  return(mldatahandler)

}

#' Run the complete meta-learning pipeline
#'
#' First we obtain the models calling to genModels and then
#' we pass those models to the metaLearning function, which will do the rest
#'
#' @param workPath String with your work path
#' @param handlerML string with the location of the handler that contains the *.dataForML files with the data
#' @param trainSpeed type of models we want to run in our metalearning function, can be FAST, FURIOUS, ALL OR BOOSTED
#' @param includeGeno bool indicating if we are going to include genotype in the meta-learning dataset
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#'
#' @return A handler with the meta learning results
#' @export
#'
#' @examples
#' handlerMeta = metaLaunch("/home/rafael", "FAST", F, "median", 30)
metaLaunch = function(workPath, trainSpeed, includeGeno, imputeMissingData, gridSearch){

  handlerMLdata = fromGenoToMLdata(workPath)

  algs = c("lda","glm", "glmnet","nb", "xgbLinear", "earth", "svmRadial",
           "rf", "bayesglm","xgbTree", "xgbDART", "C5.0Tree")

  models = NULL
  models = genModels(algs, handlerMLdata, imputeMissingData, workPath, gridSearch)

  algs = NULL
  algs[["ALL"]] = c("dnn","lda","glm", "nnet","C5.0","glmnet","nb", "xgbLinear", "earth", "svmRadial", "lasso",
                    "ridge", "evtree", "xgbTree", "rf", "bayesglm", "xgbDART")
  algs[["FAST"]] =  c("glm", "nb", "nnet", "rf", "dnn", "glmnet", "xgbTree", "xgbDART")
  algs[["FURIOUS"]] = c("glm", "nb", "nnet", "rf", "dnn", "glmnet")
  algs[["BOOSTED"]] = c("xgbLinear", "xgbTree", "xgbDART")

  handlerMeta = NULL
  handlerMeta <- metaMLtrainAndTest(models, handlerMLdata, algs, trainSpeed,imputeMissingData, workPath, includeGeno, gridSearch)
  saveRDS(handlerMeta,paste0(workPath, "/handlerMeta.rds"))

  return(handlerMeta)

}


#' Title
#'
#' @param algs Algorithms we want to use for the obtaining of the models
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param workPath String with your work path
#'
#' @return the models produced
#' @export
#'
#' @examples
#' models = genModels(c("glm","nnet","rf","xgbTree"), handler, "median", "/home/rafael", 30)
genModels = function(algs, handlerMLdata, imputeMissingData, workPath, gridSearch){

  train <- fread(handlerMLdata$train1mldata)
  train$PHENO[train$PHENO == 2] <- "DISEASE"
  train$PHENO[train$PHENO == 1] <- "CONTROL"
  ID <- train$ID
  train[,c("ID") := NULL]
  preProcValues <- preProcess(train[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  train_processed <- predict(preProcValues, train) # here we make the preprocessed values

  CVfolds <- 5
  CVrepeats <- 3
  indexPreds <- createMultiFolds(train_processed$PHENO, CVfolds, CVrepeats)
  ctrl <- trainControl(method = "repeatedcv",
                       repeats = CVrepeats,
                       number = CVfolds,
                       returnResamp = "all",
                       savePredictions = "all",
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       index = indexPreds)

  models = NULL
  for(alg in algs){
    model = NULL
    tryCatch(model <- train(PHENO ~ .,
                            data = train_processed,
                            method = alg,
                            trControl = ctrl,
                            tuneLength = gridSearch,
                            metric = "ROC")
             ,
             error = function(e){
               cat("Error when using caret algorithm",alg,":\n")
               print(e)
             })

    if(!is.null(model))
      saveRDS(model,paste0(workPath,"/models-",alg,".rds"))
    models[[alg]] = model
  }

  return(models)

}



#' Meta - Learning function
#'
#' Creating the meta machine learning data and doing
#' the learning on this new dataset.
#'
#' @param models Models obtained in the function genModels
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param algs the algorithms we are going to use to train the meta-dataset
#' @param trainSpeed type of models we want to run in our metalearning function
#' @param workPath String with your work path
#' @param includeGeno bool indicating if we are going to include genotype in the meta-learning dataset
#' @param gridSearch int with the gridSearch for the caret train method
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#'
#' @return a handler with all the meta-learning results
#' @export
#'
#' @examples
#' handlerMeta = metaMLtrainAndTest(models, handler, c("glm","nnet","rf","xgbTree"), "/home/rafael", F, 30, "median")
metaMLtrainAndTest = function(models,
                              handlerMLdata,
                              algs,
                              trainSpeed = "ALL",
                              workPath = "/home/rafa/MetaLearning/metaML/Modelos/",
                              includeGeno = F,
                              gridSearch = 30,
                              imputeMissingData = "median"){


  handlerMeta <- NULL

  checkVariantNames(handlerMLdata$train1mldata,handlerMLdata$test1mldata)
  dataTest <- fread(handlerMLdata$test1mldata)

  #same preprocess applied to the training
  dataTest$PHENO[dataTest$PHENO == 2] <- "DISEASE"
  dataTest$PHENO[dataTest$PHENO == 1] <- "CONTROL"
  ID <- dataTest$ID
  dataTest[,c("ID") := NULL]
  preProcValues <- preProcess(dataTest[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  test_processed <- predict(preProcValues, dataTest) # here we make the preprocessed values

  #dataframe which will contain the predictions from each model
  metaSet = data.frame(ID, test_processed$PHENO)
  colnames(metaSet) <- c("ID", "PHENO")

  for (i in 1:length(models)){
    preds <- predict(models[[i]], newdata = test_processed)
    colName <- paste0("Pred",i)
    metaSet[colName] <- preds
  }

  if(includeGeno){
    genoIdx <- tail(which(startsWith(colnames(test_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, test_processed[,genoIdx:length(test_processed)])
  }

  handlerMeta$metaSet = metaSet
  metaSet$ID = NULL

  #separate the new dataframe into train and test
  trainIdx<- createDataPartition(metaSet[[c("PHENO")]],p=0.75, list = FALSE, times = 1)
  trainMD<-metaSet[trainIdx,]
  testMD<-metaSet[-trainIdx,]

  trainControl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    classProbs = TRUE,
    summaryFunction = twoClassSummary)

  modelsRan <- NULL
  set.seed(1234)

  for (model in algs[[trainSpeed]]){

    modi <- train(PHENO ~ ., data = trainMD,
                  method = model,
                  trControl = trainControl,
                  tuneLength = gridSearch,
                  metric = "ROC")

    modelsRan[[model]] <- modi

  }

  ### now compare best models
  sink(file = paste(workPath, "metaML-methodPerformance.tab",sep =""), type = c("output"))
  methodComparisons <- resamples(modelsRan)
  print(summary(methodComparisons))
  sink()
  sink(file = paste(workPath, "metaML-methodTimings.tab",sep =""), type = c("output"))
  print(methodComparisons$timings)
  sink()

  handlerMeta$modelsRan = modelsRan
  handlerMeta$methodComparisons = resamples(modelsRan)

  ## pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
  ROCs <- as.matrix(methodComparisons, metric = methodComparisons$metric[1])
  meanROCs <- as.data.frame(colMeans(ROCs, na.rm = T))
  meanROCs$method <- rownames(meanROCs)
  names(meanROCs)[1] <- "meanROC"
  bestFromROC <- subset(meanROCs, meanROC == max(meanROCs$meanROC))
  bestAlgorithm <- paste(bestFromROC[1,2])
  write.table(bestAlgorithm, file = paste(workPath, "metaML-bestModel.algorithm",sep = ""), quote = F, row.names = F, col.names = F) # exports "method" option for the best algorithm
  bestModel <- modelsRan[[bestAlgorithm]]

  handlerMeta$bestAlgorithm = bestAlgorithm
  handlerMeta$bestModel = bestModel

  metaResults <- NULL
  metaResults$PHENO <- testMD$PHENO
  metaResults$predicted <- predict(bestModel, testMD)
  metaResults$probDisease <- predict(bestModel, testMD, type = "prob")[2]
  metaResults$diseaseBinomial <- ifelse(testMD$PHENO == "DISEASE", 1, 0)
  metaResults$predictedBinomial <- ifelse(metaResults$predicted == "DISEASE", 1, 0)
  write.table(metaResults, file = paste(workPath, "metaMLPredictions.tab",sep =""), quote = F, sep = "\t", row.names = F)
  handlerMeta$metaPredictions = metaResults

  confMat <- confusionMatrix(data = as.factor(metaResults$predicted), reference = as.factor(metaResults$PHENO), positive = "DISEASE")
  handlerMeta$confMat = confMat


  return(handlerMeta)
}



checkVariantNames = function(traindata,testdata){

  cat("Checking possible change of main allele between train ",traindata,
      " and test data ",testdata,"\n")
  train = fread(traindata, header = T)
  test = fread(testdata,header=T)
  snps = grep(":",colnames(train))
  snps = colnames(train)[snps]
  sentinel = F
  for(snp in snps){
    if(!(snp %in% colnames(test))){
      cat(snp," from traning data is not at test data\n")
      thesplit = stringr::str_split(snp,"_")
      commonpart = thesplit[[1]][[1]]
      letter = thesplit[[1]][[2]]
      testpos = grep(commonpart,colnames(test))
      wrongsnp = colnames(test)[testpos]
      test[,testpos] = 2 - test[,..testpos]

      cat("Converting ",wrongsnp," into ",paste0(commonpart,"_",letter),"\n")
      colnames(test)[testpos] = paste0(commonpart,"_",letter)
      sentinel = T
    }
  }
  if(sentinel){
    cat("Saving ",testdata," again to disk\n")
    fwrite(test, file = testdata, quote = F, sep = "\t", row.names = F, na = NA)
  }
  else
    cat("No need to change any variant name, identical column nanes for test/train data\n")

}


