

#' Generation of multiple repositories, given only one set of geno files
#'
#' @param workPath String with your work path
#' @param path2Geno Path to your genotype data
#' @param path2Covs Path to your covariates data
#' @param path2plink Path to PLINK
#' @param nrepos Number of repositories you want to crate
#'
#' @return
#' @export
#'
#' @examples genPartitionFromGeno("/home/rafa/",
#'                                "/home/users/gsit/juanbot/JUAN_SpanishGWAS/UNRELATED.SPAIN4.HARDCALLS.Rsq0.8",
#'                                "/home/users/gsit/juanbot/JUAN_SpanishGWAS/COVS_SPAIN",
#'                                 nrepos=3)
genPartitionFromGeno = function(workPath,
                                path2Geno = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/UNRELATED.SPAIN4.HARDCALLS.Rsq0.8",
                                path2Covs = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/COVS_SPAIN",
                                path2plink = "",
                                nrepos){

  lines <- system(paste0("wc -l ", path2Covs, ".cov | awk '{print $1}'"), intern = TRUE)
  lines = strtoi(lines)
  linesFile = floor(lines/nrepos) +1
  command = paste0("split --numeric=1 -l ",linesFile," -a1 --additional-suffix=.cov ", path2Covs,".cov ",workPath,"/covsRepo")
  mySystem(command)

  for (i in 1:nrepos){
    # inserting the first line of the covs file in every subdivided covs file
    if (i!=1) {
      firstLine = system(paste0("head -n1 ",path2Covs,".cov"), intern=TRUE)
      command = paste0("sed -i '1i",firstLine,"' ",workPath,"/covsRepo",i,".cov")
      mySystem(command)
    }

    command = paste0("awk '{ print $1,$2 }' ",workPath,"/covsRepo",i,".cov > ",workPath,"/idsRepo",i,".txt")
    mySystem(command)

    command = paste0("plink --bfile ", path2Geno, " --keep ",workPath,"/idsRepo",i,".txt --make-bed --out ",workPath,"/genoRepo",i)
    mySystem(command)

    command = paste0("awk '{ print $1,$2,$5 }' ",workPath,"/covsRepo",i,".cov > ",workPath,"/phenoRepo",i,".pheno")
    mySystem(command)

  }
}




#' Function that obtains, from the genotype data, a dataset with the data
#' prepared for the Machine Learning part, after doing the variable selection
#'
#' @param workPath String with your work path
#' @param path2Geno Path to your genotype data
#' @param path2Covs Path to your covariates data
#' @param predictor What you want to predict
#' @param path2GWAS Path to your GWAS
#' @param path2PRSICE Path to PRSICE
#' @param path2plink Path to PLINK
#' @param path2pheno Path to phenotype data (optional)
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
                            path2PRSice = "/home/users/gsit/juanbot/genoml-core/otherPackages/",
                            path2plink = "",
                            path2Pheno = paste0(workPath,"/MyPhenotype")
                            ){

  #Our first handler will always start by holding the genotype file, the covariates, the id and familial id columns
  #at the covariates, the variable name from which to predict how we want the phenotype file to be started
  h = getHandlerToGenotypeData(geno=path2Geno,
                               covs=path2Covs,
                               id="IID",
                               fid="FID",
                               predictor=predictor,
                               #With this we assure everything will be written under workPath
                               pheno=path2Pheno)

  # #Generate a holdout partition. The training part will be used to generate three ML models
  # #The test part will be used to evaluate the models (the handler is saved for later, not used in this script)
  holdout = getPartitionsFromHandler(genoHandler=h,
                                     workPath = workPath,
                                     path2plink=path2plink,
                                     how="holdout",
                                     p=0.50)

  holdout = genDataFromHandler(holdout,lazy=T)
  saveRDS(holdout,paste0(workPath,"/holdout.rds"))
  holdout = readRDS(paste0(workPath,"/holdout.rds"))

  handlerSNPs = mostRelevantSNPs(handler=getHandlerFromFold(handler=holdout,type="train",index=1),
                                 path2plink=path2plink,
                                 gwas="RISK_noSpain.tab",
                                 #gwasDef=" --beta --snp MarkerName --A1 Allele1 --A2 Allele2 --stat Effect --se StdErr --pvalue P-value",
                                 path2GWAS=path2GWAS,
                                 PRSiceexe="PRSice_linux",
                                 path2PRSice=path2PRSice,
                                 clumpField = "P-value",
                                 SNPcolumnatGWAS = "MarkerName")

  saveRDS(handlerSNPs,paste0(workPath,"/handlerSNPs.rds"))
  handlerSNPs = readRDS(paste0(workPath,"/handlerSNPs.rds"))

  mldatahandler = fromSNPs2MLdata(handler=holdout,
                                  addit="NA",
                                  path2plink=path2plink,
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
#' In a first step, we generate two datasets, for train and test,
#' and then we obtain the best model applying ML to the train dataset, for each repository.
#'
#' Once we have the best models, we run the meta learning on them to generate a more reliable result
#'
#' @param workPath String with your work path
#' @param includeGeno bool indicating if we are going to include genotype in the meta-learning dataset
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param nrepos Number of repositories we are working with
#' @param ncores Number of cores you want to use for the ML-training
#'
#' @return A handler with the meta learning results
#' @export
#'
#' @examples
#' handlerMeta = metaLaunch("/home/rafael", F, "median", 30, 3, 22)
metaLaunch = function(workPath, includeGeno, imputeMissingData, gridSearch, nrepos, ncores){

  handlersML = list()
  repoModels = list()

  for(i in 1:nrepos){

    lworkPath = paste0(workPath,"/Repo_",i)
    dir.create(lworkPath)
    handlersML[[paste0("Repo",i)]] = fromGenoToMLdata(lworkPath,
                                                      path2Geno = paste0("/home/rafajorda/exps/genPartitionFromGeno/genoRepo",i),
                                                      path2Covs = paste0("/home/rafajorda/exps/genPartitionFromGeno/covsRepo",i),
                                                      predictor = "DISEASE",
                                                      path2GWAS = "/home/rafajorda/packages/",
                                                      path2PRSice = "/home/rafajorda/packages/",
                                                      path2plink = "/home/rafajorda/packages/")

    algs = c("lda","glm", "glmnet","nb", "C5.0Tree", "earth", "svmRadial",
             "rf", "bayesglm","xgbTree", "xgbDART", "xgbLinear")

    repoModels[[paste0("Repo",i)]] = genModels(algs,
                                               handlerMLdata = handlersML[[paste0("Repo",i)]],
                                               imputeMissingData,
                                               gridSearch,
                                               ncores)
  }

  handlerMeta <- metaMLtrainAndTest(repoModels, handlersML, alg = "xgbTree", workPath, includeGeno, gridSearch, imputeMissingData)
  saveRDS(handlerMeta,paste0(workPath, "/handlerMeta"))

  return(handlerMeta)

}


#' Title
#'
#' @param algs Algorithms we want to use for the obtaining of the models
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for training
#'
#' @return the models produced
#' @export
#'
#' @examples
#' models = genModels(c("glm","nnet","rf","xgbTree"), handler, "median", "/home/rafael", 30)
genModels = function(algs, handlerMLdata, imputeMissingData, gridSearch, ncores){

  trainRepo <- fread(handlerMLdata$train1mldata)
  trainRepo$PHENO[trainRepo$PHENO == 2] <- "DISEASE"
  trainRepo$PHENO[trainRepo$PHENO == 1] <- "CONTROL"
  ID <- trainRepo$ID
  trainRepo[,c("ID") := NULL]
  preProcValues <- preProcess(trainRepo[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  train_processed <- predict(preProcValues, trainRepo) # here we make the preprocessed values


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

  ### to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)

  modelsRan = NULL
  for(alg in algs){
    model = NULL
    set.seed(1234)
    tryCatch(model <- train(PHENO ~ .,
                            data = train_processed,
                            method = alg,
                            trControl = ctrl,
                            tuneLength = gridSearch,
                            metric = "ROC")
             ,
             error = function(e){
               #cat("Error when using caret algorithm",alg,":\n")
               print(e)
             })

    if(!is.null(model)){
      modelsRan[[alg]] = model
    }
  }

  ### shut down multicore
  stopCluster(cluster)
  registerDoSEQ()

  repo = list()
  methodComparisons = resamples(modelsRan)
  repo$methodComparisons = methodComparisons

  ## pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
  ROCs <- as.matrix(methodComparisons, metric = methodComparisons$metric[1])
  meanROCs <- as.data.frame(colMeans(ROCs, na.rm = T))
  meanROCs$method <- rownames(meanROCs)
  names(meanROCs)[1] <- "meanROC"
  bestFromROC <- subset(meanROCs, meanROC == max(meanROCs$meanROC))
  bestAlgorithm <- paste(bestFromROC[1,2])
  bestModel <- modelsRan[[bestAlgorithm]]

  repo$bestAlgorithm = bestAlgorithm
  repo$bestModel = bestModel
  return(repo)

}



#' Meta - Learning function
#'
#' Creating the meta machine learning data and doing
#' the learning on this new dataset.
#'
#' @param repoModels List containing best models for each repository
#' @param handlersML List containing the ML-handlers associated to each one of the repositories
#' @param alg the algorithm we are going to use to train the meta-dataset
#' @param workPath String with your work path
#' @param includeGeno bool indicating if we are going to include genotype in the meta-learning dataset
#' @param gridSearch int with the gridSearch for the caret train method
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#'
#' @return a handler with all the meta-learning results
#' @export
#'
#' @examples
#' handlerMeta = metaMLtrainAndTest(repoModels, handlersML, "xgbTree", "/home/rafael", F, 30, "median")
metaMLtrainAndTest = function(repoModels,
                              handlersML,
                              alg = "xgbTree",
                              workPath = "/home/rafa/MetaLearning/metaML/Modelos/",
                              includeGeno = F,
                              gridSearch = 30,
                              imputeMissingData = "median"){

  metaSets <- vector(mode = "list")
  testSets <- vector(mode = "list")

  for (i in 1:length(handlersML)){
    handlerRepo = handlersML[[paste0("Repo",i)]]
    checkVariantNames(handlerRepo$train1mldata,handlerRepo$test1mldata)
    dataTest <- fread(handlerRepo$test1mldata)

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
    testSets[[i]] <- test_processed

    # paste predictions into the new df
    for (j in 1:length(repoModels)){
      aux = test_processed
      diff1 = setdiff(colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData[,-1]),colnames(aux))
      if (length(diff1) > 0){
        new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
        colnames(new_cols) = diff1
        aux = cbind(aux, new_cols)
      }
      diff2 = setdiff(colnames(aux[,-2]),colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData))
      if (length(diff2) > 0){
        aux[,diff2] = NULL
      }

      preds <- predict(repoModels[[paste0("Repo",j)]]$bestModel, newdata = aux)
      colName <- paste0("Pred",j)
      metaSet[colName] <- preds
    }

    metaSets[[i]] <- metaSet
  }

  library(dplyr)
  metaSet = bind_rows(metaSets) # new df with predictions for all repositories

  # now we add the genotype
  # first we find common columns between SNPs in each repository
  common_cols <- Reduce(intersect, lapply(testSets, colnames))
  genoIdx <- tail(which(startsWith(common_cols,"PC")), n=1) +1
  common_cols = common_cols[genoIdx:length(common_cols)]

  # once we have them, we join them by rows
  test = bind_rows(testSets)
  geno_cols = test[,..common_cols]

  # our final set is formed binding the dataframe with the predictions to the one with the genotype
  metaSet = cbind(metaSet, geno_cols)

  handlerMeta = NULL
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

  model <- train(PHENO ~ ., data = trainMD,
                  method = alg,
                  trControl = trainControl,
                  tuneLength = gridSearch,
                  metric = "ROC")

  handlerMeta$model = model

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


