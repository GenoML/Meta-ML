# Meta-ML

Meta-ML is a subproject of Geno-ML <https://genoml.github.io/>. It is an attempt to outline advanced machine learning techniques for the problems addressed in the Geno-ml platform. This is the core package of Meta-ML. This repo is under development. For any question on the pipeline please contact Rafael Jordá (rafajorda.rj@gmail.com) and Juan A. Botía (juanbotiablaya@gmail.com). 

# What this package is currently for

One of the functions of this package is make the user unaware that the Geno-ML works with plink files. So it basically makes it easy to deal with sampling on genotype data files so the genotype can be prepared for Machine Learning faster and easier. 

But of course, we enable a multi-ml model approach to PRS with ML as the first step towards federated learning.
We aim to develop Meta-Learning approaches to deal with multiple resources and to mimic (improve???) the results that a simple ML model would have on a single repository.

The steps are:

1. Generate n folds of the data, with correct case/control stratification.
2. Use external feature selection outcome (main SNPs to use) to generate the ML dataset for each fold.
3. Generate a ML model for each fold, keep all the models
4. For each ML model M do
     4.1. Interrogate M with all training examples
5. Generate a meta-ML dataset that, for each individual i we have
     5.1. {(g(i),m1(i),m2(i),...,mn(i)),status}
     5.2. status is the disease status, mj(i) is the prediction of j-th model for the ith individual
     5.3. g(i) is the SNP genotype for SNPs obtained in Step 2
6. Generate a ML model on this data
7. Evaluate on the test data


## Install the development version from GitHub:

Simply do this
```r
devtools::install_github('rafajm7/Meta-ML')
```

## Generation of ML data 
This generation involves a process of variable selection on the initial data and a conversion of the genotype files to a table in R with the variables selected.

## Obtaining of the models
In this part, we apply simple machine learning to obtain models trained with the 50% of the data generated on step 1.

## Meta Learning 
When all the models are obtained, we use them to predict on our test data (the other 50%), and with those predictions we create our new meta-dataset. Genotype data can be added to this new dataset (highly recommended), which will produce a much bigger table but with better results. 

# Credits

The development of this suite of packages is leaded by Rafael Jordá and Juan A. Botía. The Meta-Learning part has been developed by Rafael Jorda. The conversion from genotype files to ML data has been done by Juan A. Botía.

And if you use the resource please cite us, with the GitHub URL.
     
