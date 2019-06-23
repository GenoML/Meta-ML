

#' Main function that applies the Meta Learning nSmE
#' and returns the final results of the evaluation in the test set
#'
#' @param workPath orkPath String with your work path
#' @param genos Vector of strings with the name of the geno files of each repo
#' @param covs Vector of strings with the name of the covs files of each repo
#' @param pgTest String with the name of the geno files for the cohort used in test
#' @param pcTest String with the name of the covs files for the cohort used in test
#' @param testFolder String with the folder in which the test results are going to be stored
#' @param predictors Vector of strings with the name of the predictor variable of each repo
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for the ML-training
#' @param nrepos Number of repositories used in the nSmE
#'
#' @return the final results on the test data 
#' @export
#'
#' @examples finalResults = metanSmE("/home/rafa/", genos = c("UNRELATED.SPAIN4.HARDCALLS.Rsq0.8","PDBP_V3", "HBS_V3"), 
#'                                   covs = c( "COVS_SPAIN", "PDBP_covariates", "HBS_covariates" ),  "PPMI_V3", "PPMI_covariates.cov",
#'                                   "median", 5, 20, 3)
metanSmE = function(workPath, genos, covs, predictors, pgTest, pcTest, testFolder, imputeMissingData, gridSearch, ncores, nrepos){
  
  handlersML = list()
  repoModels = list()


  algsML = c("glmnet", "C5.0Tree", "earth", "svmRadial", "rf", "xgbTree")

  # for each repo, we obtain its mldata and its expert
  for (i in 1:nrepos){

    lworkPath = paste0(workPath,"/Repo_",i)
    dir.create(lworkPath)
    handlersML[[paste0("Repo",i)]] = fromGenoToMLdata(lworkPath,
                                                        path2Geno = paste0("/home/rafajorda/data/",genos[i]),
                                                        path2Covs = paste0("/home/rafajorda/data/",covs[i]),
                                                        predictor = predictors[i],
                                                        path2GWAS = "/home/rafajorda/packages/",
                                                        path2PRSice = "/home/rafajorda/packages/",
                                                        path2plink = "/home/rafajorda/packages/",
                                                        snpsSpain = handlersML[["Repo1"]]$snpsToPull,
                                                        iter = i)

    repoModels[[paste0("Repo",i)]] = genModels.nSmE(lworkPath,
                                                    algsML,
                                                    handlerMLdata = handlersML[[paste0("Repo",i)]],
                                                    imputeMissingData,
                                                    gridSearch,
                                                    ncores)

    saveRDS(repoModels[[paste0("Repo",i)]], paste0(lworkPath, "/repoModel",i,".rds"))
    saveRDS(handlersML[[paste0("Repo",i)]], paste0(lworkPath, "/handlerRepo",i,".rds"))

  }

  saveRDS(repoModels, paste0(workPath, "repoModels.rds"))
  saveRDS(handlersML, paste0(workPath,"handlersML.rds"))

  
  lworkPath = paste0(workPath,"/META/")
  dir.create(lworkPath)

  algsMML = c("C5.0Tree", "rf", "xgbTree")
  modelsMeta = NULL
  # obtain the meta learning model for each subtype
  for (t in c("MGSA", "MG", "MSA",  "MP")){
    modelsMeta[[t]] <- trainAndTestMML.nSmE(repoModels, handlersML, algs = algsMML, lworkPath, gridSearch, imputeMissingData, ncores, type = t)
  }

  lworkPath = paste0(workPath,"/",testFolder,"/")
  dir.create(lworkPath)
  
  # generate mldata for the test repository
  handlerTest = prepareFinalTest(workPath =  lworkPath,
                                       path2Geno = paste0("/home/rafajorda/data/",pgTest),
                                       path2Covs = paste0("/home/rafajorda/data/",pcTest),
                                       predictor = "PHENO_PLINK",
                                       snpsToPull = handlersML[["Repo1"]]$snpsToPull)

  
  finalResults = NULL
  # obtain final evaluation results
  for (t in c("MGSA", "MG", "MSA", "MP")){
    #modelMeta <- readRDS(paste0(workPath,"/META/modelMeta-tipo",t,".rds"))
    finalResults[[t]] = finalTest.nSmE(workPath = lworkPath, modelsMeta[[t]], handlerTest, repoModels, type = t)
    saveRDS(finalResults[[t]], paste0(workPath, "/", testFolder, "/finalResults-tipo",t,".rds"))
  }

  return(finalResults)
}






#' Main function that applies the Meta Learning 1SmE
#' and returns the final results of the evaluation in the test set
#'
#'
#' @param workPath String with your work path
#' @param pgTest String with the name of the geno files for the cohort used in test
#' @param pcTest String with the name of the covs files for the cohort used in test
#' @param testFolder String with the folder in which the test results are going to be stored
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for the ML-training
#' 
#' @return the final results on the test data (PPMI)
#' @export
#'
#' @examples finalResults = metanSmE("/home/rafa/", "PPMI_V3", "PPMI_covariates.cov",
#'                                   "median", 5, 20, 3)
meta1SmE = function(workPath, pgTest, pcTest, testFolder, imputeMissingData, gridSearch, ncores){
  
  lworkPath = paste0(workPath,"/dataRepo/")
  dir.create(lworkPath)

  handlerML = fromGenoToMLdata(lworkPath,
                               path2Geno = "/home/rafajorda/data/UNRELATED.SPAIN4.HARDCALLS.Rsq0.8",
                               path2Covs = "/home/rafajorda/data/COVS_SPAIN",
                               predictor = "DISEASE",
                               path2GWAS = "/home/rafajorda/packages/",
                               path2PRSice = "/home/rafajorda/packages/",
                               path2plink = "/home/rafajorda/packages/",
                               iter = 1)


  saveRDS(handlerML, paste0(lworkPath, "/handler.rds"))

  algsML = c("glmnet", "C5.0Tree", "earth", "svmRadial", "rf", "xgbTree")

  # generate k folds and obtain their experts
  models = genModels.1SmE(lworkPath,
                          algs = algsML,
                          handlerMLdata = handlerML,
                          imputeMissingData,
                          gridSearch,
                          ncores, k = 5)

  saveRDS(models, paste0(lworkPath, "/models.rds"))

  lworkPath = paste0(workPath,"/META/")
  dir.create(lworkPath)

  algsMML = c("C5.0Tree", "rf", "xgbTree")
  modelsMeta = NULL
  
  # obtain the meta learning model for each subtype
  for (type in c("MGSA", "MG", "MSA", "MP")){
    modelsMeta[[type]] <- trainAndTestMML.1SmE(models, handlerML, algs = algsMML, lworkPath, gridSearch, imputeMissingData, ncores, type)
  }

  
  handlerML = readRDS(paste0(lworkPath, "/handler.rds"))
  models = readRDS(paste0(lworkPath, "/models.rds"))
  
  lworkPath = paste0(workPath,"/",testFolder,"/")
  dir.create(lworkPath)
  
  # generate mldata for the test repository
  handlerTest = prepareFinalTest(workPath =  lworkPath,
                                 path2Geno = paste0("/home/rafajorda/data/",pgTest),
                                 path2Covs = paste0("/home/rafajorda/data/",pcTest),
                                 predictor = "PHENO_PLINK",
                                 snpsToPull = handlerML$snpsToPull)
  
  finalResults = NULL
  # obtain final evaluation results
  for (t in c("MGSA", "MG", "MSA", "MP")){
    #modelMeta <- readRDS(paste0(workPath,"/META2/modelMeta-tipo",t,".rds"))
    finalResults[[t]] = finalTest.1SmE(workPath = lworkPath, modelsMeta[[t]], handlerTest, models, type = t)
    saveRDS(finalResults[[t]], paste0(workPath, "/", testFolder, "/finalResults-tipo",t,".rds"))
  }
  return(finalResults)
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
#' @param path2Pheno Path to phenotype data (optional)
#' @param snpsSpain SNPs extracted from the spanish repo that will be use to extract variables in the rest of repos
#' @param iter Number that distinguish if we are going to select variables from the spanish dataset (iter = 1)
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
#'                                   "/home/users/gsit/juanbot/genoml-core/otherPackages/",
#'                                   NULL, 1)
fromGenoToMLdata = function(workPath,
                            path2Geno = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/UNRELATED.SPAIN4.HARDCALLS.Rsq0.8",
                            path2Covs = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/COVS_SPAIN",
                            predictor = "DISEASE",
                            path2GWAS = "/home/users/gsit/juanbot/JUAN_SpanishGWAS/toJuanNov7th2018/",
                            path2PRSice = "/home/users/gsit/juanbot/genoml-core/otherPackages/",
                            path2plink = "",
                            path2Pheno = paste0(workPath,"/MyPhenotype"),
                            snpsSpain = NULL,
                            iter){

  #Our first handler will always start by holding the genotype file, the covariates, the id and familial id columns
  #at the covariates, the variable name from which to predict how we want the phenotype file to be started
  h = getHandlerToGenotypeData(geno=path2Geno,
                               covs=path2Covs,
                               id="IID",
                               fid="FID",
                               predictor=predictor,
                               #With this we assure everything will be written under workPath
                               pheno=path2Pheno)

  # Generate a holdout partition. The training part will be used to generate three ML models
  # The test part will be used to evaluate the models (the handler is saved for later, not used in this script)
  holdout = getPartitionsFromHandler(genoHandler=h,
                                     workPath = workPath,
                                     path2plink=path2plink,
                                     how="holdout",
                                     p=0.75)

  holdout = genDataFromHandler(holdout,lazy=T)
  saveRDS(holdout,paste0(workPath,"/holdout.rds"))
  holdout = readRDS(paste0(workPath,"/holdout.rds"))

  # if the current repo is the spanish, apply PRSice to reduce variables
  if (iter == 1){
    handlerSNPs = mostRelevantSNPs(handler=getHandlerFromFold(handler=holdout,type="train",index=1),
                                   path2plink=path2plink,
                                   gwas="RISK_noSpain_MAF0.05.tab",
                                   path2GWAS=path2GWAS,
                                   PRSiceexe="PRSice_linux",
                                   path2PRSice=path2PRSice,
                                   clumpField = "P-value",
                                   SNPcolumnatGWAS = "MarkerName")
  }
  else{ # if not, extract only the SNPs selected from the spanish repo
    handlerSNPs = holdout
    addit = "NA"
    geno = basename(handlerSNPs$geno)
    pheno = basename(handlerSNPs$pheno)
    cov = basename(handlerSNPs$covs)
    path2Genotype=paste0(dirname(handlerSNPs$geno),"/")
    prefix=paste0("g-",geno,"-p-",pheno,"-c-",cov,"-a-",addit)
    fprefix = paste0(workPath,"/",prefix)
    
    handlerSNPs$snpsToPull = snpsSpain
    command = paste0(path2plink,"plink --bfile ",path2Genotype,geno," --extract ",
                     handlerSNPs$snpsToPull," --recode A --out ",fprefix,".reduced_genos")
    mySystem(command)
    # exports SNP list for extraction in validation set
    command = paste0("cut -f 1 ",handlerSNPs$snpsToPull," > ",fprefix,".reduced_genos_snpList")
    mySystem(command)
    handlerSNPs$rgenosSnpList = paste0(fprefix,".reduced_genos_snpList")
  }
  
  

  saveRDS(handlerSNPs,paste0(workPath,"/handlerSNPs.rds"))
  handlerSNPs = readRDS(paste0(workPath,"/handlerSNPs.rds"))

  # generate mldata for the repository
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


#' Function that obtains the expert of one repository in Meta Learning nSmE
#'
#' @param workPath String with your workpath
#' @param algs Algorithms we want to use for the obtaining of the models
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for training
#'
#' @return the best model for the repository
#' @export
#'
#' @examples
#' models = genModels.nSmE("/home/rafael",c("glm","nnet","rf","xgbTree"), handler, "median", 30, 20)
genModels.nSmE = function(workPath, algs, handlerMLdata, imputeMissingData, gridSearch, ncores){

  trainRepo <- fread(handlerMLdata$train1mldata)
  trainRepo = trainRepo[trainRepo$AGE != -9,]
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
    
  # to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)
    
  # do the training with all the algorithms included in the variable algs
  models = NULL
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
      models[[alg]] = model
    }
    saveRDS(model, paste0(workPath, "/alg-",alg,".rds"))
  }
  
  # shut down multicore
  stopCluster(cluster)
  registerDoSEQ()
  
  repo = list()
  methodComparisons = resamples(models)
  repo$methodComparisons = methodComparisons

  # pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
  ROCs <- as.matrix(methodComparisons, metric = methodComparisons$metric[1])
  meanROCs <- as.data.frame(colMeans(ROCs, na.rm = T))
  meanROCs$method <- rownames(meanROCs)
  names(meanROCs)[1] <- "meanROC"
  bestFromROC <- subset(meanROCs, meanROC == max(meanROCs$meanROC))
  bestAlgorithm <- paste(bestFromROC[1,2])
  bestModel <- models[[bestAlgorithm]]
  
  dfroc <- as.data.frame(ROCs)
  IC <- t.test(dfroc[bestAlgorithm], y = NULL, conf.level = 0.95)$conf.int
  bestMean <- mean(dfroc[[bestAlgorithm]])
  
  repo$IC = IC
  repo$bestMean = bestMean

  repo$bestAlgorithm = bestAlgorithm
  repo$bestModel = bestModel
  return(repo)

}

#' Function that generates k folds from the repo used in 1SmE and obtains the best models for each fold
#'
#' @param workPath String with your workpath
#' @param algs Algorithms we want to use for the obtaining of the models
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for training
#' @param k Number of folds generated from the repo
#'
#' @return the best models produced, for each fold
#' @export
#'
#' @examples 
#' models = genModels.1SmE("/home/rafael",c("glm","nnet","rf","xgbTree"), handler, "median", 30, 20, 5)
genModels.1SmE = function(workPath, algs, handlerMLdata, imputeMissingData, gridSearch, ncores, k){
  trainRepo <- fread(handlerMLdata$train1mldata)
  trainRepo = trainRepo[trainRepo$AGE != -9,]
  trainRepo$PHENO[trainRepo$PHENO == 2] <- "DISEASE"
  trainRepo$PHENO[trainRepo$PHENO == 1] <- "CONTROL"
  ID <- trainRepo$ID
  trainRepo[,c("ID") := NULL]
  preProcValues <- preProcess(trainRepo[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  train_processed <- predict(preProcValues, trainRepo) # here we make the preprocessed values
  
  myfolds = createFolds(y=train_processed$PHENO,k=k)
  CVfolds <- 5
  CVrepeats <- 3
  ctrl <- trainControl(method = "repeatedcv",
                         repeats = CVrepeats,
                         number = CVfolds,
                         returnResamp = "all",
                         savePredictions = "all",
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary)
  
  # for each fold, do the training with all the algorithms included in the variable algs
  finalModels = NULL
  for (i in 1:k){
      # to begin tune first begin parallel parameters
      library("parallel")
      library("doParallel")
      cluster <- makeCluster(ncores) # convention to leave 1 core for OS
      registerDoParallel(cluster)
      models = NULL
      
      for(alg in algs){
        model = NULL
        set.seed(1234)
        tryCatch(model <- train(PHENO ~ .,
                                data = train_processed[myfolds[[i]],], # data of the i-th fold
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
          models[[alg]] = model
        }
        saveRDS(model, paste0(workPath, "/fold",i,"-alg-",alg,".rds"))
      }
      
      # shut down multicore
      stopCluster(cluster)
      registerDoSEQ()
      
      repo = list()
      methodComparisons = resamples(models)
      repo$methodComparisons = methodComparisons
      
      # pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
      ROCs <- as.matrix(methodComparisons, metric = methodComparisons$metric[1])
      meanROCs <- as.data.frame(colMeans(ROCs, na.rm = T))
      meanROCs$method <- rownames(meanROCs)
      names(meanROCs)[1] <- "meanROC"
      bestFromROC <- subset(meanROCs, meanROC == max(meanROCs$meanROC))
      bestAlgorithm <- paste(bestFromROC[1,2])
      bestModel <- models[[bestAlgorithm]]
      
      dfroc <- as.data.frame(ROCs)
      IC <- t.test(dfroc[bestAlgorithm], y = NULL, conf.level = 0.95)$conf.int
      bestMean <- mean(dfroc[[bestAlgorithm]])
      
      repo$IC = IC
      repo$bestMean = bestMean
      
      repo$bestAlgorithm = bestAlgorithm
      repo$bestModel = bestModel
      
      finalModels[[i]] = repo
      
  }
  
  
  return(finalModels)
}


#' Meta - Learning function nSmE
#'
#' Creating the meta machine learning data and doing
#' the learning on this new dataset.
#'
#' @param repoModels List containing best models for each repository
#' @param handlersML List containing the ML-handlers associated to each one of the repositories
#' @param algs the algorithms we are going to use to train the meta-dataset
#' @param workPath String with your work path
#' @param gridSearch int with the gridSearch for the caret train method
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param ncores Number of cores you want to use for training
#' @param type The subtype of meta learning we are working in 
#'
#' @return the model generated by training the metaset
#' @export
#'
#' @examples
#' handlerMeta = metaMLtrainAndTes.nSmE(repoModels, handlersML, c("xgbTree", "rf"), "/home/rafael", F, 30, "median", "MGSA")
trainAndTestMML.nSmE = function(repoModels,
                                handlersML,
                                algs,
                                workPath = "/home/rafa/MetaLearning/metaML/Modelos/",
                                gridSearch = 5,
                                imputeMissingData = "median",
                                ncores,
                                type = "MP"){

  metaSets <- vector(mode = "list")
  predSets <- vector(mode = "list")

  # for each repository, obtain the predictions of every expert on evaluation data and generate its metaset
  for (i in 1:length(handlersML)){
    handlerRepo = handlersML[[paste0("Repo",i)]]
    checkVariantNames(handlerRepo$train1mldata,handlerRepo$test1mldata)

    dataTrain <- as.data.frame(fread(handlerRepo$train1mldata))
    dataTest <- as.data.frame(fread(handlerRepo$test1mldata))
    
    dataRepo <- dataTest
    #same preprocess applied to the training
    #remove individuals with invalid AGE
    dataRepo = dataRepo[dataRepo$AGE != -9,]
    dataRepo$PHENO[dataRepo$PHENO == 2] <- "DISEASE"
    dataRepo$PHENO[dataRepo$PHENO == 1] <- "CONTROL"
    ID <- dataRepo$ID
    dataRepo$ID = NULL
    preProcValues <- preProcess(dataRepo[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
    repo_processed <- predict(preProcValues, dataRepo) # here we make the preprocessed values

    #dataframe which will contain the predictions from each model
    metaSet = data.frame(ID, repo_processed$PHENO)
    colnames(metaSet) <- c("ID", "PHENO")
    predSets[[i]] <- repo_processed
    
    # remove variables in the evaluation data and not in the data used for training each expert
    # add variables set to zero present in each expert but not in the evaluation data
    for (j in 1:length(repoModels)){
      aux = repo_processed
      diff1 = setdiff(colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData[,-1]),colnames(aux))
      if (length(diff1) > 0){
        new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
        colnames(new_cols) = diff1
        aux = cbind(aux, new_cols)
      }
      diff2 = setdiff(colnames(aux[,-1]),colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData))
      if (length(diff2) > 0){
        aux[,diff2] = NULL
      }
      preds <- predict(repoModels[[paste0("Repo",j)]]$bestModel, newdata = aux) # paste predictions into the metaset
      colName <- paste0("Pred",j)
      metaSet[colName] <- preds
    }
    
    metaSets[[i]] <- metaSet
  }

  metaSet = bind_rows(metaSets) # final metaset with predictions for all repositories
  # first we find common columns between SNPs in each repository
  common_cols <- Reduce(intersect, lapply(predSets, colnames))
  
  if (type == "MGSA") { # we add age, sex and geno
    common_cols = common_cols[!common_cols %in% c("PHENO")]
    common_cols = common_cols[!startsWith(common_cols,"PC")]
    # once we have them, we join them by rows
    preds = bind_rows(predSets)
    geno_cols = preds[,common_cols]
    
    # our final set is formed binding the dataframe with the predictions to the one with the genotype
    metaSet = cbind(metaSet, geno_cols)
  }
  else if (type == "MG"){ # we add geno
    genoIdx <- tail(which(startsWith(common_cols,"PC")), n=1) +1
    common_cols = common_cols[genoIdx:length(common_cols)]
    # once we have them, we join them by rows
    preds = bind_rows(predSets)
    geno_cols = preds[,common_cols]
    
    # our final set is formed binding the dataframe with the predictions to the one with the genotype
    metaSet = cbind(metaSet, geno_cols)
    
  }
  else if (type == "MSA"){ # we add age and sex
    genoIdx <- tail(which(startsWith(common_cols,"PC")), n=1) +1
    common_cols = common_cols[1:(genoIdx-1)]
    common_cols = common_cols[!startsWith(common_cols,"PC")]
    common_cols = common_cols[!common_cols %in% c("PHENO")]
    # once we have them, we join them by rows
    preds = bind_rows(predSets)
    geno_cols = preds[,common_cols]
    
    # our final set is formed binding the dataframe with the predictions to the one with the genotype
    metaSet = cbind(metaSet, geno_cols)
  }
  else print("We don't add any extra features to the metaset")
  
  ID <- metaSet$ID
  metaSet$ID = NULL
  saveRDS(metaSet, paste0(workPath, "/metaSet-tipo",type,".rds"))

  trainIdx<- createDataPartition(metaSet[[c("PHENO")]],p=0.75, list = FALSE, times = 1)
  trainMD<-metaSet[trainIdx,]
  testMD<-metaSet[-trainIdx,]
  

  CVfolds <- 5
  CVrepeats <- 3
  indexPreds <- createMultiFolds(trainMD$PHENO, CVfolds, CVrepeats)
  ctrl <- trainControl(method = "repeatedcv",
                       repeats = CVrepeats,
                       number = CVfolds,
                       returnResamp = "all",
                       savePredictions = "all",
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       index = indexPreds)

  
  
  
  # to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)

  # finally train the metaset with all the algorithms included in the variable algs
  models = NULL
  for(alg in algs){
    model = NULL
    set.seed(1234)
    tryCatch(model <- train(PHENO ~ .,
                            data = trainMD,
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
      models[[alg]] = model
    }
    saveRDS(model, paste0(workPath, "/alg-",alg,".rds"))
  }

  # shut down multicore
  stopCluster(cluster)
  registerDoSEQ()
  
  methodComparisons = resamples(models)

  # pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
  ROCs <- as.matrix(methodComparisons, metric = methodComparisons$metric[1])
  meanROCs <- as.data.frame(colMeans(ROCs, na.rm = T))
  meanROCs$method <- rownames(meanROCs)
  names(meanROCs)[1] <- "meanROC"
  bestFromROC <- subset(meanROCs, meanROC == max(meanROCs$meanROC))
  bestAlgorithm <- paste(bestFromROC[1,2])
  bestModel <- models[[bestAlgorithm]]
  
  model <- bestModel
  metaResults <- NULL
  metaResults$PHENO <- testMD$PHENO
  metaResults$predicted <- predict(model, testMD)
  metaResults$probDisease <- predict(model, testMD, type = "prob")[2]
  metaResults$diseaseBinomial <- ifelse(testMD$PHENO == "DISEASE", 1, 0)
  metaResults$predictedBinomial <- ifelse(metaResults$predicted == "DISEASE", 1, 0)
  #write.table(metaResults, file = paste(workPath, "metaMLPredictions.tab",sep =""), quote = F, sep = "\t", row.names = F)

  confMat <- confusionMatrix(data = as.factor(metaResults$predicted), reference = as.factor(metaResults$PHENO), positive = "DISEASE")
  metaResults$confMat = confMat
  metaResults$methodComparisons = methodComparisons
  metaResults$bestAlgorithm = bestAlgorithm

  metaSet$ID = ID 
  saveRDS(model, paste0(workPath,"/modelMeta-tipo",type,".rds"))
  saveRDS(metaResults, paste0(workPath,"/metaResults-tipo",type,".rds"))
  
  return(model)
}



#' Meta - Learning function 1SmE
#'
#' Creating the meta machine learning data and doing
#' the learning on this new dataset.
#'
#' @param models Experts obtained from the function genModels
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param algs the algorithms we are going to use to train the meta-dataset
#' @param workPath String with your work path
#' @param gridSearch int with the gridSearch for the caret train method
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param ncores Number of cores you want to use for training
#' @param type The subtype of meta learning we are working in 
#'
#' @return the model generated by training the metaset
#' @export
#'
#' @examples
#' handlerMeta = metaMLtrainAndTest.1SmE(repoModels, handlersML, c("xgbTree", "rf"), "/home/rafael", F, 30, "median", "MG")
trainAndTestMML.1SmE = function(models,
                                handlerMLdata,
                                algs,
                                workPath = "/home/rafa/MetaLearning/metaML/Modelos/",
                                gridSearch = 30,
                                imputeMissingData = "median",
                                ncores,
                                type = "MP"){

  checkVariantNames(handlerMLdata$train1mldata,handlerMLdata$test1mldata)
  dataTrain <- as.data.frame(fread(handlerMLdata$train1mldata))
  dataTest <- as.data.frame(fread(handlerMLdata$test1mldata))

  totalData = dataTest
  #same preprocess applied to the training
  #remove individuals with invalid AGE
  totalData = totalData[totalData$AGE != -9,]
  totalData$PHENO[totalData$PHENO == 2] <- "DISEASE"
  totalData$PHENO[totalData$PHENO == 1] <- "CONTROL"
  ID <- totalData$ID
  totalData$ID = NULL
  preProcValues <- preProcess(totalData[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  total_processed <- predict(preProcValues, totalData) # here we make the preprocessed values

  metaSet = data.frame(ID, total_processed$PHENO)
  colnames(metaSet) <- c("ID", "PHENO")

  # obtain predictions from each expert in our evaluation data
  for (i in 1:length(models)){
    preds <- predict(models[[i]]$bestModel, newdata = total_processed)
    colName <- paste0("Pred",i)
    metaSet[colName] <- preds
  }

  if (type == "MGSA") { # we add age, sex and geno
    cols = colnames(total_processed)
    cols = cols[!cols %in% c("PHENO")]
    cols = cols[!startsWith(cols,"PC")]
    all_cols = total_processed[,cols]
    metaSet = cbind(metaSet, all_cols)
  }
  else if (type == "MG"){ # we add geno
    genoIdx <- tail(which(startsWith(colnames(total_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, total_processed[,genoIdx:length(total_processed)])
  }
  else if (type == "MSA"){ # we add age and sex
    genoIdx <- tail(which(startsWith(colnames(total_processed),"PC")), n=1) +1
    common_cols = colnames(total_processed[,1:(genoIdx-1)])
    common_cols = common_cols[!startsWith(common_cols,"PC")]
    common_cols = common_cols[!common_cols %in% c("PHENO")]

    metaSet = cbind(metaSet, total_processed[,common_cols])
  }
  else print("We don't add any extra features to the metaset")

  ID <- metaSet$ID
  metaSet$ID = NULL
  saveRDS(metaSet, paste0(workPath, "/metaSet-tipo",type,".rds"))
  #separate the new dataframe into train and test
  trainIdx<- createDataPartition(metaSet[[c("PHENO")]],p=0.75, list = FALSE, times = 1)
  trainMD<-metaSet[trainIdx,]
  testMD<-metaSet[-trainIdx,]

  CVfolds <- 5
  CVrepeats <- 3
  indexPreds <- createMultiFolds(trainMD$PHENO, CVfolds, CVrepeats)
  ctrl <- trainControl(method = "repeatedcv",
                       repeats = CVrepeats,
                       number = CVfolds,
                       returnResamp = "all",
                       savePredictions = "all",
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       index = indexPreds)

  # to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)

  # finally train the metaset with all the algorithms included in the variable algs
  models = NULL
  for(alg in algs){
    model = NULL
    set.seed(1234)
    tryCatch(model <- train(PHENO ~ .,
                            data = trainMD,
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
      models[[alg]] = model
    }
    #saveRDS(model, paste0(workPath, "/alg-",alg,"-type-", type, ".rds"))
  }
  
  # shut down multicore
  stopCluster(cluster)
  registerDoSEQ()
  
  saveRDS(models, paste0(workPath, "/models-type-",type,".rds"))
  methodComparisons = resamples(models)

  # pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
  ROCs <- as.matrix(methodComparisons, metric = methodComparisons$metric[1])
  meanROCs <- as.data.frame(colMeans(ROCs, na.rm = T))
  meanROCs$method <- rownames(meanROCs)
  names(meanROCs)[1] <- "meanROC"
  bestFromROC <- subset(meanROCs, meanROC == max(meanROCs$meanROC))
  bestAlgorithm <- paste(bestFromROC[1,2])
  bestModel <- models[[bestAlgorithm]]
  
  model <- bestModel

  metaResults <- NULL
  metaResults$PHENO <- testMD$PHENO
  metaResults$predicted <- predict(model, testMD)
  metaResults$probDisease <- predict(model, testMD, type = "prob")[2]
  metaResults$diseaseBinomial <- ifelse(testMD$PHENO == "DISEASE", 1, 0)
  metaResults$predictedBinomial <- ifelse(metaResults$predicted == "DISEASE", 1, 0)
  #write.table(metaResults, file = paste(workPath, "metaMLPredictions.tab",sep =""), quote = F, sep = "\t", row.names = F)

  confMat <- confusionMatrix(data = as.factor(metaResults$predicted), reference = as.factor(metaResults$PHENO), positive = "DISEASE")
  metaResults$confMat = confMat
  metaResults$methodComparisons = methodComparisons
  metaResults$bestAlgorithm = bestAlgorithm

  metaSet$ID = ID
  saveRDS(model, paste0(workPath,"/modelMeta-tipo",type,".rds"))
  saveRDS(metaResults, paste0(workPath,"/metaResults-tipo",type,".rds"))
  return(model)



}


#' Function that generates the ML data for the repo reserved for the test part
#'
#' @param workPath String with your work path
#' @param path2Geno Path to your genotype data
#' @param path2Covs Path to your covariates data
#' @param predictor What you want to predict
#' @param addit 
#' @param snpsToPull SNPs selected from the spanish repo, to use them to extract variables in the testing repo
#' @param path2plink Path to plink
#'
#' @return
#' @export
#'
#' @examples mldatahandler = prepareFinalTest(workPath =  "/home/rafajorda/exps2/PPMI/", 
#'                                            path2Geno = "/home/rafajorda/data/PLINK_HARD.PPMI",
#'                                            path2Covs = "/home/rafajorda/data/PPMI_covariates",
#'                                            predictor = "PHENO_PLINK",
#'                                            addit = "NA", 
#'                                            snpsToPull,
#'                                            path2plink = "/home/rafajorda/packages/")
prepareFinalTest = function(workPath =  "/home/rafajorda/exps2/PPMI/", 
                            path2Geno = "/home/rafajorda/data/PLINK_HARD.PPMI",
                            path2Covs = "/home/rafajorda/data/PPMI_covariates",
                            predictor = "PHENO_PLINK",
                            addit = "NA", 
                            snpsToPull,
                            path2plink = "/home/rafajorda/packages/"){
  
  command = paste0("cp ", path2Covs, ".cov ", workPath)
  mySystem(command)
  
  handler = getHandlerToGenotypeData(geno=path2Geno,
                                     covs=path2Covs,
                                     id="IID",
                                     fid="FID",
                                     predictor = predictor,
                                     pheno=paste0(workPath,"/MyPhenotype"))
  
  geno = basename(handler$geno)
  pheno = basename(handler$pheno)
  cov = basename(handler$covs)
  path2Genotype=paste0(dirname(handler$geno),"/")
  prefix=paste0("g-",geno,"-p-",pheno,"-c-",cov,"-a-",addit)
  fprefix = paste0(workPath,"/",prefix)
  
  handler$snpsToPull = snpsToPull # set the SNPs to pull to the ones selected in the spanish set
  command = paste0(path2plink,"plink --bfile ",path2Genotype,geno," --extract ",
                   handler$snpsToPull," --recode A --out ",fprefix,".reduced_genos")
  mySystem(command)
  # exports SNP list for extraction in validation set
  command = paste0("cut -f 1 ",handler$snpsToPull," > ",fprefix,".reduced_genos_snpList")
  mySystem(command)
  handler$rgenosSnpList = paste0(fprefix,".reduced_genos_snpList")

  # we generate de .dataForML file
  mldatahandler = fromSNPs2MLdata(handler,
                                  addit="NA",
                                  path2plink)
  
  saveRDS(mldatahandler, paste0(workPath, "/mldatahandler.rds"))
  return(mldatahandler)
}


#' Function that obtains the prediction in the test repository in the Meta Learning nSmE.
#'
#' @param workPath String with your work path 
#' @param modelMeta Model returned by the MML nSmE
#' @param handlerTest Handler with the mldata of the test set
#' @param repoModels Experts used in the generation of the metaset
#' @param type Subtype of MML that will be applied to obtain the prediction
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#'
#' @return
#' @export
#'
#' @examples finalResults = finalTest.nSmE("/home/rafa/", modelMeta, handler, repoModels, "MGSA", "median")
finalTest.nSmE = function(workPath, modelMeta, handlerTest, repoModels, type, imputeMissingData = "median"){
  
  dataFinal <- as.data.frame(fread(handlerTest$train1mldata))

  #same preprocess applied to the training
  dataFinal$PHENO[dataFinal$PHENO == 2] <- "DISEASE"
  dataFinal$PHENO[dataFinal$PHENO == 1] <- "CONTROL"
  ID <- dataFinal$ID
  dataFinal$ID = NULL
  preProcValues <- preProcess(dataFinal[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  final_processed <- predict(preProcValues, dataFinal) # here we make the preprocessed values
  
  #dataframe which will contain the predictions from each model
  metaSet = data.frame(ID, final_processed$PHENO)
  colnames(metaSet) <- c("ID", "PHENO")
  
  # obtain the metaset for the test data
  for (j in 1:length(repoModels)){
    aux = final_processed
    diff1 = setdiff(colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData[,-1]),colnames(aux))
    if (length(diff1) > 0){
      new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
      colnames(new_cols) = diff1
      aux = cbind(aux, new_cols)
    }
    diff2 = setdiff(colnames(aux[,-1]),colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData))
    if (length(diff2) > 0){
      aux[,diff2] = NULL
    }
    preds <- predict(repoModels[[paste0("Repo",j)]]$bestModel, newdata = aux)
    colName <- paste0("Pred",j)
    metaSet[colName] <- preds
  }
  
  
  if (type == "MGSA") { # we add age, sex and geno
    cols = colnames(final_processed)
    cols = cols[!cols %in% c("PHENO")]
    cols = cols[!startsWith(cols,"PC")]
    all_cols = final_processed[,cols]
    metaSet = cbind(metaSet, all_cols)
  }
  else if (type == "MG"){ # we add geno
    genoIdx <- tail(which(startsWith(colnames(final_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, final_processed[,genoIdx:length(final_processed)])
  }
  else if (type == "MSA"){ # we add age and sex
    genoIdx <- tail(which(startsWith(colnames(final_processed),"PC")), n=1) +1
    common_cols = colnames(final_processed[,1:(genoIdx-1)])
    common_cols = common_cols[!startsWith(common_cols,"PC")]
    common_cols = common_cols[!common_cols %in% c("PHENO")]
    
    metaSet = cbind(metaSet, final_processed[,common_cols])
  }
  else print("We don't add any extra features to the metaset")
  
  saveRDS(metaSet, paste0(workPath, "/metaSet-ANTES-FINAL-tipo",type,".rds"))
  
  
  ID <- metaSet$ID
  metaSet$ID = NULL
  aux = metaSet
  
  # remove variables in the metaset and not in the data used for training the MML model
  # add variables set to zero present in the MML model but not in the data
  diff1 = setdiff(colnames(modelMeta$trainingData[,-1]),colnames(aux))
  if (length(diff1) > 0){
    new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
    colnames(new_cols) = diff1
    aux = cbind(aux, new_cols)
  }
  diff2 = setdiff(colnames(aux[,-1]),colnames(modelMeta$trainingData))
  if (length(diff2) > 0){
    aux[,diff2] = NULL
    
  }
  
  metaSet = aux
  saveRDS(metaSet, paste0(workPath, "/metaSetFINAL-tipo",type,".rds"))
  
  # obtain final results
  finalResults = NULL
  finalResults$PHENO <- metaSet$PHENO
  finalResults$predFinal <- predict(modelMeta, metaSet)
  finalResults$probDisease <- predict(modelMeta, metaSet, type = "prob")[2]
  finalResults$diseaseBinomial <- ifelse(metaSet$PHENO == "DISEASE", 1, 0)
  finalResults$predictedBinomial <- ifelse(finalResults$predFinal == "DISEASE", 1, 0)
  
  confMat <- confusionMatrix(data = as.factor(finalResults$predFinal), reference = as.factor(finalResults$PHENO), positive = "DISEASE")
  finalResults$confMat = confMat
  
  return(finalResults)
}



#' Function that obtains the prediction in the test repository in the Meta Learning 1SmE.
#'
#' @param workPath String with your work path 
#' @param modelMeta Model returned by the MML nSmE
#' @param handlerTest Handler with the mldata of the test set
#' @param models Experts used in the generation of the metaset
#' @param type Subtype of MML that will be applied to obtain the prediction
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#'
#' @return
#' @export
#'
#' @examples finalResults = finalTest.nSmE("/home/rafa/", modelMeta, handler, models, "MGSA", "median")
finalTest.1SmE = function(workPath, modelMeta, handlerTest, models, type, imputeMissingData = "median"){
  
  dataFinal <- as.data.frame(fread(handlerTest$train1mldata))
  
  #same preprocess applied to the training
  dataFinal$PHENO[dataFinal$PHENO == 2] <- "DISEASE"
  dataFinal$PHENO[dataFinal$PHENO == 1] <- "CONTROL"
  ID <- dataFinal$ID
  dataFinal$ID = NULL
  preProcValues <- preProcess(dataFinal[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  final_processed <- predict(preProcValues, dataFinal) # here we make the preprocessed values
  
  #dataframe which will contain the predictions from each model
  metaSet = data.frame(ID, final_processed$PHENO)
  colnames(metaSet) <- c("ID", "PHENO")
 

  aux = final_processed
  # remove variables in the test data and not in the data used for training each expert 
  # add variables set to zero present in each expert but not in the test data
  diff1 = setdiff(colnames(models[[1]]$bestModel$trainingData[,-1]),colnames(aux)) # each expert has the same variables, so we use the first one, for example
  if (length(diff1) > 0){
    new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
    colnames(new_cols) = diff1
    aux = cbind(aux, new_cols)
  }
  diff2 = setdiff(colnames(aux[,-1]),colnames(models[[1]]$bestModel$trainingData))
  if (length(diff2) > 0){
    aux[,diff2] = NULL
  }
  for (i in 1:length(models)){
    preds <- predict(models[[i]]$bestModel, newdata = aux)
    colName <- paste0("Pred",i)
    metaSet[colName] <- preds
  }
 
  if (type == "MGSA") { # we add age, sex and geno
    cols = colnames(final_processed)
    cols = cols[!cols %in% c("PHENO")]
    cols = cols[!startsWith(cols,"PC")]
    all_cols = final_processed[,cols]
    metaSet = cbind(metaSet, all_cols)
  }
  else if (type == "MG"){ # we add geno
    genoIdx <- tail(which(startsWith(colnames(final_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, final_processed[,genoIdx:length(final_processed)])
  }
  else if (type == "MSA"){ # we add age and sex
    genoIdx <- tail(which(startsWith(colnames(final_processed),"PC")), n=1) +1
    common_cols = colnames(final_processed[,1:(genoIdx-1)])
    common_cols = common_cols[!startsWith(common_cols,"PC")]
    common_cols = common_cols[!common_cols %in% c("PHENO")]
    
    metaSet = cbind(metaSet, final_processed[,common_cols])
  }
  else print("We don't add any extra features to the metaset")
  
  saveRDS(metaSet, paste0(workPath, "/metaSet-ANTES-FINAL-tipo",type,".rds"))

  ID <- metaSet$ID
  metaSet$ID = NULL
  aux = metaSet
  
  # remove variables in the metaset and not in the data used for training the MML model
  # add variables set to zero present in the MML model but not in the data
  diff1 = setdiff(colnames(modelMeta$trainingData[,-1]),colnames(aux))
  if (length(diff1) > 0){
    new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
    colnames(new_cols) = diff1
    aux = cbind(aux, new_cols)
  }
  diff2 = setdiff(colnames(aux[,-1]),colnames(modelMeta$trainingData))
  if (length(diff2) > 0){
    aux[,diff2] = NULL
  }
  
  metaSet = aux
  saveRDS(metaSet, paste0(workPath, "/metaSetFINAL-tipo",type,".rds"))
  
  # obtain final results
  finalResults = NULL
  finalResults$PHENO <- metaSet$PHENO
  finalResults$predFinal <- predict(modelMeta, metaSet)
  finalResults$probDisease <- predict(modelMeta, metaSet, type = "prob")[2]
  finalResults$diseaseBinomial <- ifelse(metaSet$PHENO == "DISEASE", 1, 0)
  finalResults$predictedBinomial <- ifelse(finalResults$predFinal == "DISEASE", 1, 0)
  
  confMat <- confusionMatrix(data = as.factor(finalResults$predFinal), reference = as.factor(finalResults$PHENO), positive = "DISEASE")
  finalResults$confMat = confMat
  
  return(finalResults)
}



