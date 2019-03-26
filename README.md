# Meta-ML

Meta-ML is a subproject of Geno-ML <https://genoml.github.io/>. It is an attempt to outline advanced machine learning techniques for the problems addressed in the Geno-ml platform. This is the core package of Meta-ML. This repo is under development. For any question on the pipeline please contact Rafael Jordá and Juan A. Botía (juanbotiablaya@gmail.com). 

# What this package is currently for

One of the functions of this package is make the user unaware that the Geno-ML works with plink files. So it basically makes it easy to deal with sampling on genotype data files so the genotype can be prepared for Machine Learning faster and easier. 

But of course, we enable a multi-ml model approach to PRS with ML as the first step towards federated learning.
We aim to develop Meta-Learning approaches to deal with multiple resources and to mimic (improve???) the results that a simple ML model would have on a single repository.

The steps are:

1. Generate n folds of the data, with correct case/control stratification.


## Install the development version from GitHub:

Simply do this
```r
devtools::install_github('rafajm7/Meta-ML')
```

And that will be all. More help to come soon.
In the meantime, you can access the tutorials in the package.

## Generation of ML data 
This generation involves a process of variable selection on the initial data and a conversion of the genotype files to a table in R with the variables selected.

## Obtaining of the models
In this part, we apply simple machine learning to obtain models trained with the 50% of the data generated on step 1.

## Meta Learning 
When all the models are obtained, we use them to predict on our test data (the other 50%), and with those predictions we create our new meta-dataset. Genotype data can be added to this new dataset (highly recommended), which will produce a much bigger table but with better results. 

# Credits

The development of this suite of packages is leaded by Rafael Jordá and Juan A. Botía. The Meta-Learning part has been developed by Rafael Jorda. The conversion from genotype files to ML data has been done by Juan A. Botía.

And if you use the resource please cite us, with the GitHub URL.
     
