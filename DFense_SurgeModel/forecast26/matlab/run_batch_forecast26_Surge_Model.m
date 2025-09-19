% Runs Forecast 26 with Surge Model predictor for all Brazilian states
clear all; close all; clc

addpath ..\..\..\Aggregated_Data

UFv={'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA','MG','MS','MT','PA','PB','PE','PI','PR','RJ',...
    'RN','RO','RR','RS','SC','SE','SP','TO'};

for n=1:length(UFv)
    UF=UFv{n};
    forecast26_SModel_batch_def;
end
