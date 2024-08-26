This repository will allow you to generate a probability distribution with 95% confidence intervals for the number of deaths, from a set of individual risk of death predictions. It contains
1. MonteCarloAnalysisFunction.m, the actual code.
2. Mortality.csv, a file containing an example set of RODs in the correct format for input.
3. Mortality_RESULTS.fig, the output that the mortality.csv file should give when the analysis is run.
4. Monte Carlo simulation V3.docx, which describes the statistical method. 

This code works well for analysing files containing a few hundred individual RODs. More than this may take a long time. 
To analyse a bigger series of risk of deaths, there is a different version that uses multiple cores to speed up the analysis. This is available on reques
