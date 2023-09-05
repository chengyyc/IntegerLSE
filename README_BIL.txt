Notes on BayesianIntegerLogistic.R
Yingyu Cheng 09/06/2021

1. Summary
The code consists of three sections: (i) functions of Hassibi's method (with LLL) algorithms and other basic functions we will use; (ii) functions of 9 methods (3 benchmarks and 6 proposed methods) and a prediction funtion; (iii) a simple illustrative example using simulation.

To apply the code into real data, just simply run all the funtions in section (i) and (ii) first, then conduct estimation and predcition using the real data.


2. More about section (ii)

2.1. functions of 9 methods
There are 9 functions with each representing one method. The first 3 are benchmarks and last 6 are our proposed methods.
Each funtion requires three parameters: x, y and prior. Here, x is predictors, y is response, prior is the prior distribution of the coefficients which is fixed as diffused normal prior.
Each funtion will return a list consisting of coefficients, estimated_probability and classification_result.

2.2. prediction function: pred_BIL
This function requires two parameters: obj and x_test. Here, x_test is the testing set, and obj is any output of the functions above, which is a list.
This function will return the prediction result of the testing set.


3. Explaination of the functions of 9 methods in section (ii)

3.1. Benchmarks
• BILOpt: first obtain the MAP estimate through optimizing algorithm then directly round it.
• BILBayMedian: first obtain the median of posterior samples, then directly round it.
• BILBayMean: first obtain the mean of posterior samples, then directly round it.

3.2. Proposed methods
First obtain the posterior samples, with each column representing a coefficient, each row representing a posterior sample.
• BILRoundMedian, BILRoundMode, BILRoundModeRow: project all the posterior samples to integer space by rounding, then separately
- first get the median of each column, then combine them as the estimate.
- first get the mode of each column, then combine them as the estimate.
- treat each row as a whole and get the mode of the rows.
• BILLLLMedian, BILLLLMode, BILLLLModeRow: project all the posterior samples to integer space by the LLL searching algorithm, then separately
- first get the median of each column, then combine them as the estimate.
- first get the mode of each column, then combine them as the estimate.
- treat each row as a whole and get the mode of the rows.


4. Instructions on simulation

4.1. Parameters
There are many parameters in simulation, and you can always change them to the values you want.
• q: Number of predictors (not including intercept).
• qb: Number of binary predictors (if applicable).
• N: Sample size.
• Rho: Correlation between predictors.
• x0: Intercept.
• my_prior: The prior distribution of coefficients.
• B: Predictors. You can use function continuous_predictors to generate continuous predictors, use function binary_predictors to generate binary predictors, or use function con_bin_predictors to generate both continuous and binary predictors.
Within function continuous_predictors and con_bin_predictors, each continuous predictor is now fixed as generating from normal distribution N(0,1). However, you can change the argument sigma_b inside the function to change the standard deviation of the distribution of continuous predictors.
• z: Coefficients. You can generate z using any distribution you want.
• n: Number of repeatings. Since we will calculte the average results of these n repeatings, the larger n is, the more general our results are, and also, the longer it will take to run the codes.

4.2. Simulation process
There are three parts in simulation process: (i) generating the data; (ii) conducting in-sample predcition; (iii) conduction out-of-sample prediction.
The metrics used to do evaluation are: (i) MSE of coefficients; (ii) in-sample AUC; (iii) out-of-sample AUC. Since there are 9 methods in total, the calculating and recording process are kind of repetitive.





