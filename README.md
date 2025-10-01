# deg-pipeline
Matthew Notes
functions implemented:
1. rewritten new deg pipeline: streamlined for better readability (functionality still lacking for sure)
2. filter design: filters duplicate columns, overspecified columns, linear combinations within data, and 
near zero variance (all to clean the data better)
3. wrote some tests to test if these all work on very basic summarized experiment datasets. 

More information
Raw RNA data is passed into the pipeline. The goal of the pipeline is to first filter out all of the data
that is useless, and then analyze it to find anomalies that need to be looked at. RNA data often has 
cooked data in the sense that some samples will have the same data, and some traits will have little 
to no variance. it is hard to test these. The pipeline filters out these and only focuses on statistically
significant data

The desired output of most of these functions are Summarized Experiment objects. It is a container for experiment 
data where rows are features, columns are samples, with coordinated metadata. 
assays: one or more matrices with identical dimensions and dimnames
rowData: per-feature data for samples: genes, IDs, etc. Think features.
colData: per-sample metadata. think of this as the information for each sample. 
metadata: experiment level information
generally, rows are featues and columns are samples. 

A lot of functions will also create a design matrix. This is similar to a Summarized Experiment (SE) object. They 
use the data and essentialy just turn everything into numbers: rows are samples, and columns are predictors. each
column is of a different feature (group, age, etc). 

The flow for the pipeline is:
raw data gets passed in, then is pseudobulked (each unique individual in the experiment has multiple samples, we 
group them together. or we group samples by type like brain region, etc. just think of this as a big grouping
algorithm). After the pseudobulk, the data is filtered as needed. then, DEG methods for RNA-seq analysis are done,
and the results are written to excel for scientists to use. 

Ben wants the pipeline to handle errors better, in the sense where it will give the user feedback of what to change
in the data. I implemented a on and off switch between automatic and manual filtering, where manual gives feedback 
and doesn't change anything while automatic does the opposite. He also wants the writing to excel to be better, where
the written data is more organized, maybe color coded, etc. 

My advice:
Many of the methods used in this pipeline are from the R packages edgeR and deSEQ2. You could read up on documentation
if you want, but honestly I find i† prety boring and chat can do most of the explaining. I would definitely play around
with R, make some example datasets using chatGPT or cursor templates to figure out the syntax. Definitely learn the 
syntax and how the data.frame object works, and near the end play around with writexl (package for excel in R). 
I use R in vscode/cursor, but R studio is also good for writing R. i dont like it though because it drains my battery. 
I'm going to upload a screenshot of the slide I made for a presentation about this pipeline, hopefully it gives you more
clarity. 
Also talk to Ben for advice, but write stuff down and do research on what he says afterwards. He's lowkey really smart 
and knows a bunch about this R and data analysis stuff, so he's gonna drop a lot of good wisdom in the meetings but also
a lot of technical terms.