# D-FENSE_SurgeModel

1) D-FENSE (Dynamics for Epidemic Surveillance and Evaluation)

Paulo Antonio Andrade Esquef (National Laboratory for Scientific Computing)

Add other team members 


2) Repository Structure


Aggregated_Data: raw data (CSV spreadsheets) of dengue cases per state used model estimation. 

DFense_SurgeModel:  dengue cases prediction based on na average surge model 

  |_ validation1:  material related to validation 1 challenge

      |_ matlab: Matlab scritps needed to run (run_batch_v1_predictor_Surge_Model.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases

      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.


  |_ validation2:  material related to validation 2 challenge

      |_ matlab: Matlab scritps needed to run (run_batch_v2_predictor_Surge_Model.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases

      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.

  |_ validation3:  material related to validation 3 challenge

      |_ matlab: Matlab scritps needed to run (run_batch_v3_predictor_Surge_Model.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases
  
      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.

|_ forecast26:  material related dengue forecast cases for season 2026

      |_ matlab: Matlab scritps needed to run (run_batch_forecast26_Surge_Model.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases
  
      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.



3) Libraries and Dependencies


For DFense_SurgeModel: 

matlab functions: readtable.m, buffer.m (Signal Processing Toolbox), lsqcurvefit.m (Optimization Toolbox)
 

4) Data and Variables

Only the time-series of raw number of dengue cases per state along epidemic weeks have been used. Data are
available from:
 
https://github.com/americocunhajr/D-FENSE/tree/main/DengueSprint2025_DataAggregated


5) Model Training

DFense_SurgeModel: for each state (UF), time-series of raw dengue cases, in the defined range for each validation, have been organized in blocks of 52 samples (one year). Assuming that the dengue surges happen about the same time (around EW 15) each year, an average or typical surge curve has been obtained. Assuming the surge is symmetrical with respect to its local maximum, a centralized (to its peak) version of the surge is obtained. From the centralized typical surge, we estimate the parameters (L,k,x0) of the derivative of the logistic model, by means of a nonlinear estimator (lsqcurvefit.m, with algoritm 'trust-region-reflective'). Then, we use a template matching filter scheme to find the local maxima of the cross-correlation coefficient sequence, between the model surge (template) and the observed surges over time. After time-synchronizing the model with a given observed surge, we calculate the amplitude gain that, when applied to the model, matches it with the each observed surge. We do that for each surge and obtain a set of amplitude gains, which are positive. The dengue cases prediction is simply given by a gain that multiplies the surge model.  
Assuming that the set of gains follows a log-normal distribution, we use the set of gains to estimate the related mean and sigma of a log-normal distribution. To predict the dengue cases, we generate 10000 gains from the previously estimated log-normal distribution and apply it to the model surge, properly placed in time. From the set of these 10000 case predictions, the median, lower- and upper-bounds of the 50%, 80%, 0%, 90%, and 95% prediction intervals are calculated. Finally, we cropped out the predictions to be in the range from EW 41 of a given year to EW 40 of the subsequent year.     
    
     
6) References

None.

7) Data Usage Restriction

None. 

DFense_SurgeModel: from the trained/estimated typical surge model (see section 5 above): After time-synchronizing the surge model with a given observed surge, we calculate the amplitude gain that, when applied to the model, matches it with the each observed surge. We do that for each surge and obtain a set of amplitude gains, which are positive. The dengue cases prediction is simply given by a gain that multiplies the surge model. Assuming that the set of gains follows a log-normal distribution, we use the set of gains to estimate the related mean and sigma of a log-normal distribution. To predict the dengue cases, we generate 10000 gains from the previously estimated log-normal distribution and apply it to the model surge, properly placed in time. From the set of these 10000 case predictions, the median, lower- and upper-bounds of the 50%, 80%, 0%, 90%, and 95% prediction intervals are calculated. Finally, we cropped out the predictions to be in the range from EW 41 of a given year to EW 40 of the subsequent year.    


8) Predictive Uncertainty


From the set of 10000 case predictions (for each state and each validation), we used matlab function prctile.m (percentiles of a sample) to obtain the median, as well as the lower- and upper bounds of 50%, 80%, 90%, and 95% prediction intervals. The median of the case predictions is the 50% percentile. The lower-bounds for the 50%, 80%, 90%, and 95% prediction intervals are, respectively, the 25%, 10%, 5%, and 2.5% percentiles. The upper-bounds for the 50%, 80%, 90%, and 95% prediction intervals are, respectively, the 75%, 90%, 95%, and 97.5% percentiles.

Summary:
median prediction: 50% percentile
50% prediction interval: from 25% percentile to 75% percentile  
80% prediction interval: from 10% percentile to 90% percentile
90% prediction interval: from 5% percentile to 95% percentile
95% prediction interval: from 2.5% percentile to 97.5% percentile

