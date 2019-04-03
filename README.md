# Meta-ML (MML for short)

Meta-ML is a subproject of Geno-ML <https://genoml.github.io/>. It is an attempt to outline advanced machine learning techniques for the problems addressed in the Geno-ml platform. This is the core package of Meta-ML. This repo is under development. For any question on the pipeline please contact Rafael Jordá (rafajorda.rj@gmail.com) and Juan A. Botía (juanbotiablaya@gmail.com).

Meta machine-learning (MML) is an approach to ML in which both the ML method we use (search bias) and the model hypotesis we optimize (representational bias) are not fixed and can be changed along the training process. One of the types of MML we use here is inductive MML. In this approach, a pool of machine learning algorithms is used on the learning data. We combine learning data with predictions from each model (i.e. the experts and their opinions on samples) to generate a MML dataset so we transform a MML problem into a ML problem. In this manner, we incorporate many search heuristics as information on how they perform on the ML problem and the machine can learn about it to improve the prediction.

Basically, if our learning data is {(xi,yi): 1 <= i <= n} and the xi refers to the genotype data (0, 1, 2) for the i-th SNP we consider in the learning process, inductive MML generates a new dataset {(xi,ei,yi): 1 <= i <= n} where ei=(ei1,ei2,...,ein) such that eij is the opinion of expert j-th on the i-th example.



# Install the development version from GitHub:

Simply do this
```r
devtools::install_github('rafajm7/Meta-ML')
```

# Inductive MML


We aim to develop MML approaches to deal with multiple resources and to mimic (improve???) the results that a simple ML model would have on a single repository. What we outline here as an algorithm is the MML, 1 source, many experts. We can have following types of inductive MML:

* __1SmE__: the type just mentioned. There is only one dataset, many possible experts working on folds for the same dataset, a MML dataset with opinion from all experts on the dataset and and a final model from that dataset.

* __nSmE__: there are many datasets, many possible experts. For each dataset, we use Caret (or any other tool avaiable) to select the best model. Thus, each dataset generates a best expert. All the original data is used to make predictions with all experts. Thus, a N samples dataset generates a MML dataset with N samples and M additional predictors, each for the j-th expert prediction to the corresponding sample.

Both these algorithms can have version with genotype as predictors in the META-ML dataset. They would be written as __1SmEg__ and __nSmEg__ ML2 (meta-ml) algorithms.

## An algorithm for 1SmE

In this section we offer a sketch of an inductive meta-learning algorithm to be applied on a single dataset. 
The steps for the __1SmE__ algorithm are

```r
0. Consider the problem of generating a model for predicting a polygenic risk score on a case/control set, 
  from plink files of cases and controls. The following steps correspond to a meta-learning approach for 
  the generation of the prediction model.
1. Separate data into 75% for training (Tr), 25% for evaluation (Ev) with case/control stratification.
2. Do feature selection on the Tr dataset (i.e. select the main SNPs to use). Once SNPs are selected, 
  filter the columns of the training data keeping only those SNPs,  to generate the training dataset.
3. Let Alg={alg1, alg2, ..., algn} the different ML techniques
  (e.g. neural networks, decision trees) we want to use for the generation of basic models.
  Generate appropriate folds on Tr for each basic ML model.
4. For each ML algorithm in Alg, and its training data chunk do

     4.1. Use Caret + the algorithm + training data to generate
      the best possible model m

5. For each ML model M do

     5.1. Interrogate M with all training examples

6. Generate a Meta-ML (MML) dataset such that, for each subject i in the Tr set,  we generate a
  new data example

     6.1. {(g(i),m1(i),m2(i),...,mn(i),status)} where status is the disease
      status for the subject, mj(i) is the prediction of j-th model
      for the i-th subject. And g(i) is the SNP genotype for SNPs obtained in
      Step 2 for the i-th individual
      

7. Let MetaAlg the technique we want to use for the meta-learner. We have then
  to use Caret + MetaAlg + MML data from point 6.1., to generate the best
  possible model m
8. Evaluate on the test data
```

First ideas on these kinds of approaches were tested on the StatLog project
in the eighties, check this report (Feng,C., Sutherland,A., King,S.,
  Muggleton,S. & Henery,R. (1993). Comparison of Machine Learning Classifiers
  to Statistics and Neural Networks. AI & Stats Conf. 93. )  
  <https://doi.org/10.1080/08839519508945477>. We´re now revisiting
  these ideas: more data availability and more powerful heuristics can deliver
  a different outcome now.

## An algorithm for nSmE

In this section we offer the specification of an algorithm in which there are n repositories of genetic data. Possibly, these individuals will have different genetic backgrounds across repositories, they will vary in size and will have different biological covariates.

Consider the problem of generating a model for predicting a polygenic risk score on a number of different repositories. Because of the particularities of the problem, we have to keep the repositories separated.  The following steps correspond to a meta-learning approach for the generation of a prediction model capable of maintaining a good error while learning particularities from all repositories segregated. As repositories are separated, the algorithm will have a concurrent part and a centralised part. We start with the concurrent part and finish gathering all results from the different runs into a single model.

The steps for the __nSmE__ algorithm are

```r
0. Let n be the number of repositories (plink files + covariates + phenotype) we have
1. Do this in parallel: 
  for r in {1..n} do
  1.1. Separate data of the r-th repository, Dr, into 75% for training (Tr,r), 
    25% for evaluation (Ev,r) with case/control stratification.
  1.2. Do feature selection on the Tr,r dataset (i.e. select the main SNPs to use). 
    Once SNPs are selected, filter the columns of the training data keeping 
    only those SNPs,  to generate the training dataset.
  1.3. Let Alg={alg1, alg2, ..., algn} a set of ML techniques
  (e.g. neural networks, decision trees) we use for the generation of basic models.   
  1.4. Use Caret + Alg + Tr,r to generate the best possible model m.
  1.5. Let Mml_r be the ml model for the r-th repository. Let also SNPml_r be the set 
    of SNPs that were selected as the main SNPs from the genotype to be used in the 
    model generation process. This parallel run returns (SNPml_r, Mml_r)

2. Let {(SNPml_i, Mml_i): 1 <= i <= n} be the results gathered from all repositories. 
  We have now to create a meta-learning dataset to learn from all repositories. 
  2.1 For each model m in the set {(SNPml_i, Mml_i)} do
    2.1.1. For each repository r in {1..n} do
      Interrogate m with all individuals in Ev,r
      Interrogate m with all individuals in Tr,r
      Let m(g(i,r)) be the prediction of the model m on the genotype of the i-th individual of repository r.
    end 2.1.1
  end 2.1
  2.2. Let us define 
    MMLD={(m1(g(i,r)), m2(g(i,r)), m3(g(i,r)), m4(g(i,r)), p1(i), p2(i),..., pk(i),status(i)), 1 <= i <= s} 
    to be the meta-learning dataset in which the mj(g(i,r)) elements refer to the prediction of the 
    model created from repository j, on individual i which belongs ro repository r. 
    The p´s are additional covariates on the individual (e.g. age, gender, genetic PCAs).

3. Use Caret + Alg + MMLD to generate the final predictor
```

There are three main steps: (1) parallel repository model learning, (2) meta-learning dataset creation and (3) final model learning on meta-learning data. 


# How to use the API to implement MML

One of the features this package has is trying to make it transparent to the user that we are using plink files. So it basically makes it easy to deal with creating appropriate data samples on genotype data files so the Machine Learning experiments can be designed faster and easier.

This data generation step involves a process of variable selection on the initial data and a conversion of the genotype files to a tibble in R with the variables selected. With the code provided we can do it. We´ll illustrate it with an example.

## Working with the plink files without working with the plink files

Let us suppose we want to work with the bed/bim/fam plink files that we can find at a specific folder, `~/mymldata/`. We asume there are only a set of genotype files there but we can also specify data per file. We will use the following algorithms (they are caret algorithm names) as basic learners: lda, glm, glmnet, nb, xgbLinear, earth, svmRadial, rf, bayesglm, xgbTree, xgbDART, C5.0Tree. Have a look at
<https://topepo.github.io/caret/available-models.html> for all available techniques.

Let us suppose we do not want to include the genotype data into the Meta-ML data set, so we are implementing __1SmE__ and not __1SmEg__. We are going to impute data with the mediam and generate grid sizes for hyper-parameter optimisation of size 30. We can generate the basic models with this simple code.

```r
imputeMissingData = "mediam"
workPath = "~/mymldata/"
gridSearch = 30

handlerMLdata = fromGenoToMLdata(workPath)
algs = c("lda","glm", "glmnet","nb", "xgbLinear", "earth", "svmRadial",
           "rf", "bayesglm","xgbTree", "xgbDART", "C5.0Tree")
basicModels = genModels(algs, handlerMLdata, imputeMissingData, workPath, gridSearch)
```

It starts from the whole genotype, divides it into 50% for train and 50% for test, does variant selection on the 50% training data, and using these variants, it generates a ML learning dataset including covariates, population structure PCAs and the phenotype as an outcome.

More in detail, this is how we can do that.

Firstly we get a handler to the genotype data. A hander is basically an abstraction of all files we manage in the background, but it is actually a file name and samples holder to work on creating sampling strategies without having to call plink each time. Plink commands are only used when the genotype data is really required. So when we call the function, no call to plink made.

Our first handler will always start by holding the genotype file, the covariates, the id and familial id columns at the covariates, the variable name from which to predict how we want the phenotype file to be called as it will be automatically generated by `getHandlerToGenotypeData()`.


```r
h = getHandlerToGenotypeData(geno=path2Geno,
                               covs=path2Covs,
                               id="IID",
                               fid="FID",
                               predictor=predictor,
                               #With this we assure everything will be written under workPath
                               pheno=paste0(workPath,"/MyPhenotype"))
```

Then we partition the data into 50%, 50% as we said. And this generates another handler, we store at `holdout`. Again, the genotype plink files remains untouched yet.

```r
holdout = getPartitionsFromHandler(genoHandler=h,
                                     workPath = workPath,
                                     path2plink="",
                                     how="holdout",
                                     p=0.50)
```

## Selecting the relevant SNPs for Learning

Now, the `genDataFromHandler()` call really generates the new plink data files. This slightly modifies the handler so it returns a new one, that we store at the same variable.  The `mostRelevantSNPs()` uses PRISice software to select the variants based on the GWAS data and p-value thresholds.

```r
holdout = genDataFromHandler(holdout,lazy=T)

handlerSNPs = mostRelevantSNPs(handler=getHandlerFromFold(handler=holdout,type="train",index=1),
                                 path2plink="",
                                 gwas="RISK_noSpain.tab",
                                 gwasDef=" --beta --snp MarkerName --A1 Allele1 --A2 Allele2 --stat Effect --se StdErr --pvalue P-value",
                                 path2GWAS=path2GWAS,
                                 PRSiceexe="PRSice_linux",
                                 path2PRSice=path2PRSice,
                                 clumpField = "P-value",
                                 SNPcolumnatGWAS = "MarkerName")
```

Now, saving the hander will allow to remember this particular feature selection results (note that PRSice is CPU intensive so better to avoid run it again in the same data if we can; but also we save that for the sake of reproducibility).

## Generation of the ML data for the basic learners

Finally, we can generate the machine learning dataset that can be readily used with any learning algorithm. Note how we tell the `fromSNPs2MLdata()` method that it has to use the variants we selected with the calls from above, through the parameter `fsHandler`.

```r
mldatahandler = fromSNPs2MLdata(handler=holdout,
                                  addit="NA",
                                  path2plink="",
                                  fsHandler=handlerSNPs)
cat("Your train dataset is at",mldatahandler$train1mldata,"\n",
"Your test dataset is at",mldatahandler$test1mldata,"\n")
```

Once we have the training data, we can call `genModels()`. What is does is step 4 of the algorithm, as follows: it prepares the `PHENO` data column to be used as a factor, it gets rid of the ID so it is not used in the learning phase, and imputes missing data if any.

```r
train = fread(handlerMLdata$train1mldata)
train$PHENO[train$PHENO == 2] = "DISEASE"
train$PHENO[train$PHENO == 1] = "CONTROL"
ID = train$ID
train[,c("ID") := NULL]
preProcValues = preProcess(train[,-1],
  method = c(paste(imputeMissingData,"Impute", sep = "")))
train_processed = predict(preProcValues, train)
```

## Generation of the base models

At its current implementation, `genModels()` creates 5 folds on the data, and three repeats. All algorithms in Alg are tested with exactly the same data and the same repeats.

```r
CVfolds = 5
CVrepeats = 3
indexPreds = createMultiFolds(train_processed$PHENO, CVfolds, CVrepeats)
ctrl = trainControl(method = "repeatedcv",
                       repeats = CVrepeats,
                       number = CVfolds,
                       returnResamp = "all",
                       savePredictions = "all",
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       index = indexPreds)
```

As you see above, we call `createMultiFolds()` indexes on the data so we can use all algorithms in Alg with exactly same samples. Therefore, we also need only one `trainControl()` object.

Below, as you can see we call `caret::train()` method with each algorithm and keep all results.  

```r
models = NULL
for(alg in algs){
  model = NULL
  model = train(PHENO ~ .,
                          data = train_processed,
                          method = alg,
                          trControl = ctrl,
                          tuneLength = gridSearch,
                          metric = "ROC")
  saveRDS(model,paste0(workPath,"/models-",alg,".rds"))
  models[[alg]] = model
}
```
## Generation of the META-learning datasets

Once we have generated a bunch of models with the 50% of the original data, we will have the other 50% to generate the MML learning data.

First of all, we have to adapt the SNPs on the test data with the SNPs that we used on the train data. As we are separating one big cohort into two smaller ones, probably the MAF of some variants will change. We have to account for that with the method `checkVariantNames()`. Then we prepare data and impute as above.


```r
checkVariantNames(handlerMLdata$train1mldata,handlerMLdata$test1mldata)
dataTest = fread(handlerMLdata$test1mldata)

dataTest$PHENO[dataTest$PHENO == 2] = "DISEASE"
dataTest$PHENO[dataTest$PHENO == 1] = "CONTROL"
ID = dataTest$ID
dataTest[,c("ID") := NULL]
preProcValues = preProcess(dataTest[,-1],
method = c(paste(imputeMissingData,"Impute", sep = "")))
test_processed = predict(preProcValues, dataTest)
```

Now we create a new data frame, `metaSet`

```r
metaSet = data.frame(ID, test_processed$PHENO)
colnames(metaSet) = c("ID", "PHENO")
```
adding, for all individuals, a column with the prediction of all models we stored in `models` in the former code.

```r
for (i in 1:length(models)){
  preds = predict(models[[i]], newdata = test_processed)
  colName = paste0("Pred",i)
  metaSet[colName] = preds
}
```

Optionally, we include the genotype in the data (i.e. we add columns for all the SNPs we use for learning).

```r
if(includeGeno){
  genoIdx = tail(which(startsWith(colnames(test_processed),"PC")), n=1) +1
  metaSet = cbind(metaSet, test_processed[,genoIdx:length(test_processed)])
}
```
And what follows now is a normal Caret based learning process. We again use all models but this time, we make them to compete between them in figuring out the final outcome, from the outcomes of each basic model. We create a data partition, the `trainControl` object and finally we call `caret::train()`

```r
handlerMeta$metaSet = metaSet
metaSet$ID = NULL

trainIdx= createDataPartition(metaSet[[c("PHENO")]],p=0.75, list = FALSE, times = 1)
trainMD=metaSet[trainIdx,]
testMD=metaSet[-trainIdx,]

trainControl = trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

modelsRan = NULL
set.seed(1234)

for (model in algs[[trainSpeed]]){
  modelsRan[[model]] = train(PHENO ~ ., data = trainMD,
                method = model,
                trControl = trainControl,
                tuneLength = gridSearch,
                metric = "ROC")
}
```
As we generated all models with the same `trainControl` object, we can evaluate all in the same way. We get the model with best ROC. All results will be stored in the  `handlerMeta` list.

```r
methodComparisons = resamples(modelsRan)
handlerMeta$modelsRan = modelsRan
handlerMeta$methodComparisons = resamples(modelsRan)

ROCs = as.matrix(methodComparisons, metric = methodComparisons$metric[1])
meanROCs = as.data.frame(colMeans(ROCs, na.rm = T))
meanROCs$method = rownames(meanROCs)
names(meanROCs)[1] = "meanROC"
bestFromROC = subset(meanROCs, meanROC == max(meanROCs$meanROC))
bestAlgorithm = paste(bestFromROC[1,2])
bestModel = modelsRan[[bestAlgorithm]]

handlerMeta$bestAlgorithm = bestAlgorithm
handlerMeta$bestModel = bestModel
```

And finally, we get the predictions and the confusion matrix.

```r


metaResults = NULL
metaResults$PHENO = testMD$PHENO
metaResults$predicted = predict(bestModel, testMD)
metaResults$probDisease = predict(bestModel, testMD, type = "prob")[2]
metaResults$diseaseBinomial = ifelse(testMD$PHENO == "DISEASE", 1, 0)
metaResults$predictedBinomial = ifelse(metaResults$predicted == "DISEASE", 1, 0)
handlerMeta$metaPredictions = metaResults

confMat = confusionMatrix(data = as.factor(metaResults$predicted),
  reference = as.factor(metaResults$PHENO), positive = "DISEASE")
handlerMeta$confMat = confMat

```






# Credits

Part of this code has been adapted from the GenoML scripts for basic model generation. They were developed by Mike Nalls and Faraz Fahgri. The development of this suite of packages is leaded by Rafael Jordá and Juan A. Botía. The Meta-Learning part has been developed by Rafael Jorda. The conversion from genotype files to ML data has been done by Juan A. Botía.

And if you use the resource please cite us, with the GitHub URL.
