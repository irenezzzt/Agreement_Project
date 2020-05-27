
Background: 

DNA methylation can be measured with several platforms, each with its own inherent strengths, limitations and biases. To illustrate and better understand the biases associated with two platforms, we analysed and compared the estimated methylation levels from the Illumina 450K array and targeted custom capture bisulfite sequencing (TCCBS) on 42 samples from 20 individuals.

In this work, Illumina 450K data were normalized with SWAN and funtooNorm to remove probe type and colour channel-related biases. For TCCBS data, we built zero inflated negative binomial models for the methylated and unmethylated counts separately, as a function of GC content. Measurement of bias and its corresponding correction method are proposed based on the model estimated counts. We then model the discrepancy between estimated methylation levels from two platforms and ranked the importance of predictors by a random forest model.


The R scripts are for DNA methylation data on TCCBS platform. Please run in order:

1. data_prep_small.R: Clean and extract a small subset (n=5000). Run principle component analysis (PCA) on GC contents and prepare data in different formats.

2. random_forest_models.R: Model the methylation level and count by random forest models.

3. other_models.R: Model the methylation level and count by zero inflated negative binomial models with different specifications.

4. bias_correction.R: Define a ``meaningful'' prediction from the ZINB models. Obtain the bias corrected count and methylation level.

