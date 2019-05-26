

#' Main function that applies the meta-learning nSmE
#' and returns the final results of the evaluation in the test set
#'
#' @param workPath orkPath String with your work path
#' @param genos Vector of strings with the name of the geno files of each repo
#' @param covs Vector of strings with the name of the covs files of each repo
#' @param predictors Vector of strings with the name of the predictor variable of each repo
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for the ML-training
#' @param nrepos Number of repositories used in the nSmE
#'
#' @return the final results on the test data (PPMI)
#' @export
#'
#' @examples
metanSmE = function(workPath, genos, covs, predictors, imputeMissingData, gridSearch, ncores, nrepos){
  
  handlersML = list()
  repoModels = list()

  lworkPath = paste0(workPath,"/dataRepo/")
  dir.create(lworkPath)

  algsML = c("glmnet", "C5.0Tree", "earth", "svmRadial", "rf", "xgbTree")

  for (i in 1:nrepos){

    lworkPath = paste0(workPath,"/Repo_",i)
    dir.create(lworkPath)
    print(handlersML[["Repo1"]]$snpsToPull)
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
  
  #repoModels = readRDS(paste0(workPath, "/repoModels.rds"))
  #handlersML = readRDS(paste0(workPath,"/handlersML.rds"))
  
  lworkPath = paste0(workPath,"/META/")
  dir.create(lworkPath)
  
  algsMML = c("C5.0Tree", "rf", "xgbTree")

  for (t in c("ALL", "GENO", "AGE-SEX",  "PREDS")){
    print(t)
    model <- trainAndTestMML.nSmE(repoModels, handlersML, algs = algsMML, lworkPath, gridSearch, imputeMissingData, ncores, type = t)
  }

  lworkPath = paste0(workPath,"/PPMI/") #por como lo hizo juan, estoy obligado a poner el fichero PPMI_covariates.cov en esa carpeta ya antes
  
  handlerTest = prepareFinalTest(workPath =  lworkPath,
                                       path2Geno = "/home/rafajorda/data/PPMI_V3",
                                       path2Covs = "/home/rafajorda/data/PPMI_covariates",
                                       predictor = "PHENO_PLINK",
                                       snpsToPull = handlersML[["Repo1"]]$snpsToPull)

  finalResults = NULL
  for (t in c("ALL", "GENO", "AGE-SEX", "PREDS")){
    print(t)
    modelMeta <- readRDS(paste0(workPath,"/META/modelMeta-tipo",t,".rds"))
    finalResults[[t]] = finalTest(workPath = lworkPath, modelMeta, handlerTest, repoModels, type = t)
    saveRDS(finalResults, paste0(workPath, "/PPMI/finalResults-tipo",t,".rds"))
  }
  
  return(finalResults)
}






#' Main function that applies the meta-learning 1SmE
#' and returns the final results of the evaluation in the test set
#'
#'
#' @param workPath String with your work path
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for the ML-training
#' 
#' @return the final results on the test data (PPMI)
#' @export
#'
#' @examples
#' handlerMeta = metaLaunch("/home/rafael", F, "median", 30, 22)
meta1SmE = function(workPath, imputeMissingData, gridSearch, ncores){
  
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
   
  models = genModels.1SmE(lworkPath,
                          algs = algsML,
                          handlerMLdata = handlerML,
                          imputeMissingData,
                          gridSearch,
                          ncores)
   
  saveRDS(models, paste0(lworkPath, "/models.rds"))
   
  #models = readRDS("~/exps1SmE/dataRepo/models.rds")
  lworkPath = paste0(workPath,"/META/")
  dir.create(lworkPath)
  
  algsMML = c("C5.0Tree", "rf", "xgbTree")
  
  for (type in c("ALL", "GENO", "AGE-SEX", "PREDS")){
    print(type)
    model <- trainAndTestMML.1SmE(models, handlerML, algs = algsMML, lworkPath, gridSearch, imputeMissingData, ncores, type)
  }
  
  lworkPath = paste0(workPath,"/PPMI/")

  handlerTest = prepareFinalTest(workPath =  lworkPath,
                                 path2Geno = "/home/rafajorda/data/PPMI_V3",
                                 path2Covs = "/home/rafajorda/data/PPMI_covariates",
                                 predictor = "PHENO_PLINK",
                                 snpsToPull = handlerML$snpsToPull)
  
  finalResults = NULL
  for (t in c("ALL", "GENO", "AGE-SEX", "PREDS")){
    modelMeta <- readRDS(paste0(workPath,"/META/modelMeta-tipo",t,".rds"))
    finalResults[[t]] = finalTest.1SmE(workPath = lworkPath, modelMeta, handlerTest, models, type = t)
    saveRDS(finalResults, paste0(workPath,"/PPMI/finalResults-tipo",t,".rds"))
  }
  
  return(finalResults)
}



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
  
  lines <- system(paste0("wc -l ", path2Covs,".cov | awk '{print $1}'"), intern = TRUE)
  lines = strtoi(lines)
  linesFile = floor(lines/nrepos) +1

  COVS = as.data.frame(fread(paste0(path2Covs,".cov")))

  for (i in 1:nrepos){
    if (i == nrepos) sample = COVS
    else{
     sample = COVS[sample(nrow(COVS), linesFile),]
     COVS = COVS[!COVS$IID %in% sample$IID,]
    }
    write.table(sample, file = paste0(workPath,"covsRepo",i,".cov"), quote = F, sep = "\t", row.names = F)
    
    idsRepo = sample[,1:2]
    write.table(idsRepo, file = paste0(workPath,"idsRepo",i,".txt"), quote = F, sep = "\t", row.names = F)
    
    command = paste0("plink --bfile ", path2Geno, " --keep ",workPath,"/idsRepo",i,".txt --make-bed --out ",workPath,"/genoRepo",i)
    mySystem(command)
    
    phenoRepo = sample[,c(1,2,5)]
    write.table(phenoRepo, file = paste0(workPath,"/phenoRepo",i,".pheno"), quote = F, sep = "\t", row.names = F)
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
#' @param iter
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

  # #Generate a holdout partition. The training part will be used to generate three ML models
  # #The test part will be used to evaluate the models (the handler is saved for later, not used in this script)
  holdout = getPartitionsFromHandler(genoHandler=h,
                                     workPath = workPath,
                                     path2plink=path2plink,
                                     how="holdout",
                                     p=0.75)

  holdout = genDataFromHandler(holdout,lazy=T)
  saveRDS(holdout,paste0(workPath,"/holdout.rds"))
  holdout = readRDS(paste0(workPath,"/holdout.rds"))

  if (iter == 1){
    handlerSNPs = mostRelevantSNPs(handler=getHandlerFromFold(handler=holdout,type="train",index=1),
                                   path2plink=path2plink,
                                   gwas="RISK_noSpain.tab",
                                   #gwasDef=" --beta --snp MarkerName --A1 Allele1 --A2 Allele2 --stat Effect --se StdErr --pvalue P-value",
                                   path2GWAS=path2GWAS,
                                   PRSiceexe="PRSice_linux",
                                   path2PRSice=path2PRSice,
                                   clumpField = "P-value",
                                   SNPcolumnatGWAS = "MarkerName")
  }
  else{
    handlerSNPs = holdout
    addit = "NA"
    geno = basename(handlerSNPs$geno)
    pheno = basename(handlerSNPs$pheno)
    cov = basename(handlerSNPs$covs)
    path2Genotype=paste0(dirname(handlerSNPs$geno),"/")
    ### options passed from list on draftCommandOptions.txt
    prefix=paste0("g-",geno,"-p-",pheno,"-c-",cov,"-a-",addit)
    fprefix = paste0(workPath,"/",prefix)
    
    handlerSNPs$snpsToPull = snpsSpain
    command = paste0(path2plink,"plink --bfile ",path2Genotype,geno," --extract ",
                     handlerSNPs$snpsToPull," --recode A --out ",fprefix,".reduced_genos")
    mySystem(command)
    # exports SNP list for extraction in validation set
    command = paste0("cut -f 1 ",handlerSNPs$snpsToPull," > ",fprefix,".reduced_genos_snpList")
    #cat("The command",command,"\n")
    mySystem(command)
    handlerSNPs$rgenosSnpList = paste0(fprefix,".reduced_genos_snpList")
  }
  
  

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


#' Title
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
    
  ### to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)
    
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
  
  ### shut down multicore
  stopCluster(cluster)
  registerDoSEQ()
  
  repo = list()
  methodComparisons = resamples(models)
  repo$methodComparisons = methodComparisons

  ## pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
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

#' Title
#'
#' @param workPath String with your workpath
#' @param algs Algorithms we want to use for the obtaining of the models
#' @param handlerMLdata handler that contains the *.dataForML files with the data
#' @param imputeMissingData string with the type of preprocess we want to apply to the data
#' @param gridSearch int with the gridSearch for the caret train method
#' @param ncores Number of cores you want to use for training
#'
#' @return the best models produced, for each repo
#' @export
#'
#' @examples 
#' models = genModels.1SmE("/home/rafael",c("glm","nnet","rf","xgbTree"), handler, "median", 30, 20)
genModels.1SmE = function(workPath, algs, handlerMLdata, imputeMissingData, gridSearch, ncores){
  trainRepo <- fread(handlerMLdata$train1mldata)
  trainRepo = trainRepo[trainRepo$AGE != -9,]
  trainRepo$PHENO[trainRepo$PHENO == 2] <- "DISEASE"
  trainRepo$PHENO[trainRepo$PHENO == 1] <- "CONTROL"
  ID <- trainRepo$ID
  trainRepo[,c("ID") := NULL]
  preProcValues <- preProcess(trainRepo[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  train_processed <- predict(preProcValues, trainRepo) # here we make the preprocessed values
  
  k = 5
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
  
  finalModels = NULL
  for (i in 1:k){
      ### to begin tune first begin parallel parameters
      library("parallel")
      library("doParallel")
      cluster <- makeCluster(ncores) # convention to leave 1 core for OS
      registerDoParallel(cluster)
      models = NULL
      
      for(alg in algs){
        model = NULL
        set.seed(1234)
        tryCatch(model <- train(PHENO ~ .,
                                data = train_processed[myfolds[[i]],],
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
      
      ### shut down multicore
      stopCluster(cluster)
      registerDoSEQ()
      
      repo = list()
      methodComparisons = resamples(models)
      repo$methodComparisons = methodComparisons
      
      ## pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
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
#' @param alg the algorithm we are going to use to train the meta-dataset
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
#' handlerMeta = metaMLtrainAndTes.nSmE(repoModels, handlersML, "xgbTree", "/home/rafael", F, 30, "median", "ALL")
trainAndTestMML.nSmE = function(repoModels,
                                handlersML,
                                algs,
                                workPath = "/home/rafa/MetaLearning/metaML/Modelos/",
                                gridSearch = 5,
                                imputeMissingData = "median",
                                ncores,
                                type = "PREDS"){

  metaSets <- vector(mode = "list")
  predSets <- vector(mode = "list")

  for (i in 1:length(handlersML)){
    handlerRepo = handlersML[[paste0("Repo",i)]]
    checkVariantNames(handlerRepo$train1mldata,handlerRepo$test1mldata)

    dataTrain <- as.data.frame(fread(handlerRepo$train1mldata))
    dataTest <- as.data.frame(fread(handlerRepo$test1mldata))
    
    dataRepo <- rbind(dataTrain,dataTest)
    
    #dataRepo <- dataTest
    dataRepo = dataRepo[dataRepo$AGE != -9,]
    #same preprocess applied to the training
    dataRepo$PHENO[dataRepo$PHENO == 2] <- "DISEASE"
    dataRepo$PHENO[dataRepo$PHENO == 1] <- "CONTROL"
    ID <- dataRepo$ID
    dataRepo$ID = NULL
    #dataRepo[,c("ID") := NULL]
    preProcValues <- preProcess(dataRepo[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
    repo_processed <- predict(preProcValues, dataRepo) # here we make the preprocessed values

    #dataframe which will contain the predictions from each model
    metaSet = data.frame(ID, repo_processed$PHENO)
    colnames(metaSet) <- c("ID", "PHENO")
    predSets[[i]] <- repo_processed
    
    # paste predictions into the new df
    for (j in 1:length(repoModels)){
      aux = repo_processed
      diff1 = setdiff(colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData[,-1]),colnames(aux))
      if (length(diff1) > 0){
        new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
        colnames(new_cols) = diff1
        aux = cbind(aux, new_cols)
        print("ha añadido")
      }
      diff2 = setdiff(colnames(aux[,-1]),colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData))
      if (length(diff2) > 0){
        aux[,diff2] = NULL
        print("ha quitado")
      }
      preds <- predict(repoModels[[paste0("Repo",j)]]$bestModel, newdata = aux)
      colName <- paste0("Pred",j)
      metaSet[colName] <- preds
    }
    
    metaSets[[i]] <- metaSet
  }

  metaSet = bind_rows(metaSets) # new df with predictions for all repositories
  saveRDS(metaSet, paste0(workPath, "/metaSetANTES-tipo",type,".rds"))
  # first we find common columns between SNPs in each repository
  common_cols <- Reduce(intersect, lapply(predSets, colnames))
  
  if (type == "ALL") { # we add age, sex and geno
    common_cols = common_cols[!common_cols %in% c("PHENO")]
    common_cols = common_cols[!startsWith(common_cols,"PC")]
    # once we have them, we join them by rows
    preds = bind_rows(predSets)
    geno_cols = preds[,common_cols]
    
    # our final set is formed binding the dataframe with the predictions to the one with the genotype
    metaSet = cbind(metaSet, geno_cols)
  }
  else if (type == "GENO"){ # we add geno
    genoIdx <- tail(which(startsWith(common_cols,"PC")), n=1) +1
    common_cols = common_cols[genoIdx:length(common_cols)]
    # once we have them, we join them by rows
    preds = bind_rows(predSets)
    geno_cols = preds[,common_cols]
    
    # our final set is formed binding the dataframe with the predictions to the one with the genotype
    metaSet = cbind(metaSet, geno_cols)
    
  }
  else if (type == "AGE-SEX"){ # we add age and sex
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
  
  print(length(colnames(metaSet)))
  
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

  
  
  
  ### to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)

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

  ### shut down multicore
  stopCluster(cluster)
  registerDoSEQ()
  
  methodComparisons = resamples(models)

  ## pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
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
#' @param models Models obtained in the function genModels
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
#' handlerMeta = metaMLtrainAndTest.1SmE(repoModels, handlersML, "xgbTree", "/home/rafael", F, 30, "median", "GENO")
trainAndTestMML.1SmE = function(models,
                                handlerMLdata,
                                algs,
                                workPath = "/home/rafa/MetaLearning/metaML/Modelos/",
                                gridSearch = 30,
                                imputeMissingData = "median",
                                ncores,
                                type = "PREDS"){

  checkVariantNames(handlerMLdata$train1mldata,handlerMLdata$test1mldata)
  dataTrain <- as.data.frame(fread(handlerMLdata$train1mldata))
  dataTest <- as.data.frame(fread(handlerMLdata$test1mldata))
  nSep <- nrow(dataTrain) # to divide the metaset into the training and the test part
  print(nSep)
  totalData <- rbind(dataTrain, dataTest)

  #same preprocess applied to the training
  totalData = totalData[totalData$AGE != -9,]
  totalData$PHENO[totalData$PHENO == 2] <- "DISEASE"
  totalData$PHENO[totalData$PHENO == 1] <- "CONTROL"
  ID <- totalData$ID
  totalData$ID = NULL
  #totalData[,c("ID") := NULL]
  preProcValues <- preProcess(totalData[,-1], method = c(paste(imputeMissingData,"Impute", sep = ""))) # note here we pick impute method (KNN or median),  we can also exclude near zero variance predictors and correlated predictors
  total_processed <- predict(preProcValues, totalData) # here we make the preprocessed values

  metaSet = data.frame(ID, total_processed$PHENO)
  colnames(metaSet) <- c("ID", "PHENO")

  for (i in 1:length(models)){
    print(models[[i]]$bestAlgorithm)
    preds <- predict(models[[i]]$bestModel, newdata = total_processed)
    colName <- paste0("Pred",i)
    metaSet[colName] <- preds
  }

  if (type == "ALL") { # we add age, sex and geno
    cols = colnames(total_processed)
    cols = cols[!cols %in% c("PHENO")]
    cols = cols[!startsWith(cols,"PC")]
    all_cols = total_processed[,cols]
    metaSet = cbind(metaSet, all_cols)
  }
  else if (type == "GENO"){ # we add geno
    genoIdx <- tail(which(startsWith(colnames(total_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, total_processed[,genoIdx:length(total_processed)])
  }
  else if (type == "AGE-SEX"){ # we add age and sex
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
  print(dim(metaSet))
  #separate the new dataframe into train and test
  #trainIdx<- createDataPartition(metaSet[[c("PHENO")]],p=0.75, list = FALSE, times = 1)
  trainMD<-metaSet[1:nSep,]
  testMD<-metaSet[-(1:nSep),]

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

  print("llega a entrenar")
  ### to begin tune first begin parallel parameters
  library("parallel")
  library("doParallel")
  cluster <- makeCluster(ncores) # convention to leave 1 core for OS
  registerDoParallel(cluster)

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
  
  ### shut down multicore
  stopCluster(cluster)
  registerDoSEQ()
  
  saveRDS(models, paste0(workPath, "/models-type-",type,".rds"))
  methodComparisons = resamples(models)

  ## pick best model from model compare then output plots in this case, its picked via ROC, maximizing the mean AUC across resamplings
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


#' Title
#'
#' @param workPath String with your work path
#' @param path2Geno Path to your genotype data
#' @param path2Covs Path to your covariates data
#' @param predictor What you want to predict
#' @param addit 
#' @param snpsToPull 
#' @param path2plink Path to plink
#'
#' @return
#' @export
#'
#' @examples 
prepareFinalTest = function(workPath =  "/home/rafajorda/exps2/PPMI/", 
                            path2Geno = "/home/rafajorda/data/PLINK_HARD.PPMI",
                            path2Covs = "/home/rafajorda/data/PPMI_covariates",
                            predictor = "PHENO_PLINK",
                            addit = "NA", 
                            snpsToPull,
                            path2plink = "/home/rafajorda/packages/"){
  
  handler = getHandlerToGenotypeData(geno=path2Geno,
                                     covs=path2Covs,
                                     id="IID",
                                     fid="FID",
                                     predictor = predictor,
                                     pheno=paste0(workPath,"/MyPhenotype"))
  
  # for (i in 1:length(handlersML)){
  #   
  #   command = paste0("sed -i '/^$/d' ", handlersML[[i]]$snpsToPull)  #quitamos líneas en blanco
  #   mySystem(command)
  #   
  #   if (i == 1){
  #     command = paste0("cp ", handlersML[[1]]$snpsToPull, " ", workPath, "/merged.snpsToPull2")
  #     mySystem(command)
  #   }
  #   else{
  #     command = paste0("sed -i '1d' ", handlersML[[i]]$snpsToPull) # quitamos la primera línea de los que no sean el primero
  #     mySystem(command)
  #     
  #     command = paste0("cat ",handlersML[[i]]$snpsToPull, " >> ", workPath, "/merged.snpsToPull2")
  #     mySystem(command)
  #   }
  # }
  
  geno = basename(handler$geno)
  pheno = basename(handler$pheno)
  cov = basename(handler$covs)
  path2Genotype=paste0(dirname(handler$geno),"/")
  ### options passed from list on draftCommandOptions.txt
  prefix=paste0("g-",geno,"-p-",pheno,"-c-",cov,"-a-",addit)
  fprefix = paste0(workPath,"/",prefix)
  
  handler$snpsToPull = snpsToPull
  command = paste0(path2plink,"plink --bfile ",path2Genotype,geno," --extract ",
                   handler$snpsToPull," --recode A --out ",fprefix,".reduced_genos")
  mySystem(command)
  # exports SNP list for extraction in validation set
  command = paste0("cut -f 1 ",handler$snpsToPull," > ",fprefix,".reduced_genos_snpList")
  #cat("The command",command,"\n")
  mySystem(command)
  handler$rgenosSnpList = paste0(fprefix,".reduced_genos_snpList")

  # we generate de .dataForML file
  mldatahandler = fromSNPs2MLdata(handler,
                                  addit="NA",
                                  path2plink)
  
  saveRDS(mldatahandler, paste0(workPath, "/mldatahandler.rds"))
  return(mldatahandler)
}


#' Title
#'
#' @param workPath 
#' @param modelMeta 
#' @param handlerTest 
#' @param repoModels 
#' @param type 
#' @param imputeMissingData 
#'
#' @return
#' @export
#'
#' @examples
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
  
  #ojo con los PCA, hablar con Juan
  for (j in 1:length(repoModels)){
    aux = final_processed
    diff1 = setdiff(colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData[,-1]),colnames(aux))
    if (length(diff1) > 0){
      new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
      colnames(new_cols) = diff1
      aux = cbind(aux, new_cols)
      print("ha añadido")
    }
    diff2 = setdiff(colnames(aux[,-1]),colnames(repoModels[[paste0("Repo",j)]]$bestModel$trainingData))
    if (length(diff2) > 0){
      aux[,diff2] = NULL
      print("ha quitado")
    }
    preds <- predict(repoModels[[paste0("Repo",j)]]$bestModel, newdata = aux)
    colName <- paste0("Pred",j)
    metaSet[colName] <- preds
  }
  
  
  if (type == "ALL") { # we add age, sex and geno
    cols = colnames(final_processed)
    cols = cols[!cols %in% c("PHENO")]
    cols = cols[!startsWith(cols,"PC")]
    all_cols = final_processed[,cols]
    metaSet = cbind(metaSet, all_cols)
  }
  else if (type == "GENO"){ # we add geno
    genoIdx <- tail(which(startsWith(colnames(final_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, final_processed[,genoIdx:length(final_processed)])
  }
  else if (type == "AGE-SEX"){ # we add age and sex
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
  
  diff1 = setdiff(colnames(modelMeta$trainingData[,-1]),colnames(aux))
  if (length(diff1) > 0){
    new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
    colnames(new_cols) = diff1
    aux = cbind(aux, new_cols)
    print("ha añadido luego")
  }
  diff2 = setdiff(colnames(aux[,-1]),colnames(modelMeta$trainingData))
  if (length(diff2) > 0){
    print("ha quitado luego")
    print(diff2)
    aux[,diff2] = NULL
    
  }
  
  metaSet = aux
  saveRDS(metaSet, paste0(workPath, "/metaSetFINAL-tipo",type,".rds"))
  
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



#' Title
#'
#' @param workPath 
#' @param modelMeta 
#' @param handlerTest 
#' @param models 
#' @param type 
#' @param imputeMissingData 
#'
#' @return
#' @export
#'
#' @examples
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
  diff1 = setdiff(colnames(models[[1]]$bestModel$trainingData[,-1]),colnames(aux))
  if (length(diff1) > 0){
    new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
    colnames(new_cols) = diff1
    aux = cbind(aux, new_cols)
    print("ha añadido")
  }
  diff2 = setdiff(colnames(aux[,-1]),colnames(models[[1]]$bestModel$trainingData))
  if (length(diff2) > 0){
    aux[,diff2] = NULL
    print("ha quitado")
  }
  for (i in 1:length(models)){
    preds <- predict(models[[i]]$bestModel, newdata = aux)
    colName <- paste0("Pred",i)
    metaSet[colName] <- preds
  }
 
  if (type == "ALL") { # we add age, sex and geno
    cols = colnames(final_processed)
    cols = cols[!cols %in% c("PHENO")]
    cols = cols[!startsWith(cols,"PC")]
    all_cols = final_processed[,cols]
    metaSet = cbind(metaSet, all_cols)
  }
  else if (type == "GENO"){ # we add geno
    genoIdx <- tail(which(startsWith(colnames(final_processed),"PC")), n=1) +1
    metaSet = cbind(metaSet, final_processed[,genoIdx:length(final_processed)])
  }
  else if (type == "AGE-SEX"){ # we add age and sex
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
  
  diff1 = setdiff(colnames(modelMeta$trainingData[,-1]),colnames(aux))
  if (length(diff1) > 0){
    new_cols = data.frame(matrix(0, ncol = length(diff1), nrow = nrow(aux)))
    colnames(new_cols) = diff1
    aux = cbind(aux, new_cols)
    print("ha añadido luego")
  }
  diff2 = setdiff(colnames(aux[,-1]),colnames(modelMeta$trainingData))
  if (length(diff2) > 0){
    aux[,diff2] = NULL
    print("ha quitado luego")
  }
  
  metaSet = aux
  saveRDS(metaSet, paste0(workPath, "/metaSetFINAL-tipo",type,".rds"))
  
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



