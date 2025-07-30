% Tasks of this script (for a given BR state):
% 1) Estimate a surge model (logistic model) from the average surge behavior across the
% years, by means of a nolinear estimator. We can consider the surge model
% as a typical (expected) surge in a given year.

% 2) Use the surge model as a template for a matched filter to obtain the
% maxima of correlation coefficients (between the surge model and an observed surge) above 0.8.  
% The time instants of these maxima serve to time-synchronize the surge
% model and the observed surges. 

% 3) Obtain a set of gains across the years that, when applied to the surge model, match its scale to 
% that of a given observed surge (time-synchronized with the surge model).
% The behavior of the set of gains may be informative somehow. One could check 
% how these gains behave among several states (neighboring or not). One could also check 
% whether there are strong correlations between these gains with the external factors related
% to climate variables or other dimensions.

% 4) Reconstruct the time-series based on the surge model and the set of
% gains. 

% 5) Provide a 'Validation 3' Forecast based on the model of the average surge and
% a simple statistical analysis of the set of gains.  



M = readtable(['DengueSprint2025_AggregatedData_',UF,'.csv']); % reads aggregated
% data from the chosen state, from 2010 to 2025

PP=52;  % period of 52 Epidemic Weeks (EW)

% Validation 2 
ind_EW25_2024=649+52+52; % sample index of EW 25 of 2024

dcases=M{:,2};   % time-series of the surges across the EWs 

dcases_orig=dcases;  % stores the original time-series up to EW 40 of 2025

dcases=dcases(1:ind_EW25_2024);  % crops time-series at EW 25 of 2024 (including)

% From the time-series 'dcases' of the dengue surges, we obtain an average surge curve.
% The hypothesis is that the surges have a seasonality of 52 Epidemic Week
% (EW), and that a properly scaled average surge is sufficient to represent
% an observed surge in a given year.

DC=buffer(dcases,PP); % organizes vector dcases in a matrix with PP rows
% each column of DC has cases for 52 EWs. 

typ_DC=mean(DC')'; % typical dcases waveform in 52 EWs: from EW 1 to EW 52  

[aux1,ii_max]=max(typ_DC);  % finds local maximum of typ_DC (template)

ii_max=min(ii_max,25); % limits ii_max to 25, otherwise circ_lag non-positive, which cannot happen


% Before estimating the surge model (logistic), we will center the surge around its maximum

circ_lag=26-ii_max;  % circle-shifts the curve to centralize it around its maximum
typ_DC_centered=circshift(typ_DC,circ_lag);  % typical surge centered around its maximum

% Surge Model Estimation 

% Initial parameters of the logistic model
L=120000;  
k=0.3;
n0=26;  % time-shift
P=[L,k,n0];
n=0:length(typ_DC_centered)-1; n=n(:);  % time basis

% Surge Model Estimation (nonlinear) 
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
fun = @(P,n) (P(1).*P(2).*exp(P(2).*(n-(P(3)))))./(1+exp(P(2).*(n-(P(3)))).^2); % surge model
P_est = lsqcurvefit(fun,P,n,typ_DC_centered,[500;0.15;20],[370000;0.5;28],options); % find estimated model parameters 


Estimated_Model_surge=fun(P_est,n);  % Synthesize surge from estimated model surge 


% Note: here, we are estimating just one model for the average surge, which
% we suppposed to be enough to represent all other surges, up to a scale factor.   
% Thus, for the purpose of surge forecasting, we only need a surge model and some statistics
% related to the set of scale gains. 
% We could also estimate several surge models, one for each year. This way,
% we could also carry out a surge forecast using statistics
% for the surge model parameters (L, k, and n0).

% Now, we match filter the time-series dcases with the surge model (as a
% template), to obtain the cross-correlation coefficients (ccc) over time.
% Maxima of the ccc above a given threshold indicate that there is a match
% between an observed surge and the surge model template, centered at the maximum of the template. 

surge_model_lag=Estimated_Model_surge;  % template surge model

% Now, we consider 'surge_model_lag' as a template for the typical behavior
% of EW 41 (given year) to EW 40 of the subsequent year. Then, we obtain a set of gains to match
% the real surge and the template

dcases_ext_out=[dcases; zeros(27,1)]; % extends dcases with 27 zeros to facilite 
% matching of the last surge

Nt=length(surge_model_lag); % length of the surge model (52 EWs)

y1=filter(surge_model_lag,1,dcases_ext_out);  % matched filters dcases through the surge model

stn=norm(surge_model_lag,2);  % norm-2 of the surge model

y2=filter(ones(Nt,1),1,abs(dcases_ext_out).^2); % filters the square of dcases through a moving average filters.
% this is to compute the correlation coefficient (see below)

rr=y1./(stn.*sqrt(y2));  % correlation coefficients (between -1 < rr < 1)

offset=find(surge_model_lag==max(surge_model_lag));  % time index of the local maximum of the surge model

thr_rr=0.7; % threshold for the correlation coefficient 
inddd=(rr>=thr_rr);  % indices of rr above thr_rr   
rr(~inddd)=0;   % sets to zero rr < thr_rr  

% find peaks of the correlation indices above thr_rr

idx_max = find( (rr(3:end) < rr(2:end-1)) & (rr(1:end-2) < (rr(2:end-1) ) ) ) +1;
idx_max = idx_max(:);

% removes idx_max indices that are too close to each other (multiple local peaks 
% of the correlation coefficient clustered together);

dd=diff(idx_max);
ind1_dd=find(dd<20);  % finds indices of idx_max which less than 20 samples appart
if ~isempty(ind1_dd)
    for kk=1:length(ind1_dd)
        [v1,ind1]=min(rr(idx_max(ind1_dd(kk):ind1_dd(kk)+1)));
        aux(kk)=ind1_dd(kk)-1+ind1;
    end
    idx_max(aux)=[]; clear aux ind1
end

indpeak=idx_max-offset+1;  % time indices of local maxima of the correlation coeficcients  

cce=[zeros(offset,1);dcases_ext_out]; % extends dcases with offset initial zeros, to compensate
% for the time-shift applyed to the model surge template

% Calculates the surge model gains to match the amplitude of the typical
% surge model with that of the observed surges 

lip=length(indpeak); % number of local maxima in rr
g=zeros(lip,1);  % initializes with zeros vector to store the set of gains
x=surge_model_lag;   % attributes surge_model_lag to x (for notation clarity)  
for kk=1:lip
    a=cce(indpeak(kk)+(1:Nt));   % observed surge centered around its maximum
    g(kk)=a'*x./(x'*x);   % calculates amplitude gain for each observed surge
end


% Now, we can reconstruct the time-series of surges using a model of the average surge
% and the set of gains

dcases_rec=zeros(size(cce)); % initializes the time-series with zeros

% Constructs the modeled time-series from the set of gains and the typical
% surge model
for kk=1:lip
    dcases_rec(indpeak(kk)+(1:Nt))=dcases_rec(indpeak(kk)+(1:Nt))+g(kk)*surge_model_lag;
end

dcases_rec(1:offset)=[];  % removes extra initial 'offset' values of dcases_rec
dcases_rec(end-(circ_lag-1):end)=[];  % crops reconstructed cases
% up to EW 39 of 2024 (including)


% Note: dcases_rec goes to EW 39 of 2024, i.e., it includes a
% 'deterministic' forecast of EW 26 to EW 39 of 2024
forecast_determ=dcases_rec(ind_EW25_2024+1:end);


% A simple predictor to forecast the cases in Validation 3 is to generate
% a set of values of g for the next season based on the mean and variance
% of the previous values of g. For simplicity a log normal distribution can
% be assumed for generating g.

[param]=lognfit(g);

mg=param(1);  % mean of value of g (from 2010 to 2024)
sigma=param(2); % standard deviation of value of g (from 2010 to 2024)

MC=10000;  % number of Monte Carlo runs for forecast surges in 2025
% Note that forecast surge from EW 26 to 39 of 2024 may be best taken from the fitted
% surge around EW 15 of 2024.

g_MC=lognrnd(mg,sigma,MC,1);  % set of randomly generated gains from log-normal distribution


surge_model=[surge_model_lag(1:circ_lag); circshift(surge_model_lag,-circ_lag)];  
% surge model from EW 40 of 2024 to EW 52 of 2025

forecast_cases_v1=zeros(length(surge_model),MC);  % matrix to store realizations of the forecast surges 
% from EW 40 of 2024 to EW 52 of 2025

for kk=1:MC
    forecast_cases_v1(:,kk)=g_MC(kk)*surge_model;
end




% Obtain approximations of the 50%, 80%, 90%, and 95% prediction intervals (data driven)

set_prctile=[2.5 5 10 25 50 75 90 95 97.5]; % 2.5 to 97.5% percentiles
PP=prctile(forecast_cases_v1',set_prctile); % calculates the percentiles and stores in PP
forecast_cases_q2p5 = PP(1,:)';  % 2.5% percentile 
forecast_cases_q5 = PP(2,:)';  % 5% percentile
forecast_cases_q10 = PP(3,:)';  % 10% percentile
forecast_cases_q25 = PP(4,:)';  % 25% percentile
mean_cases_forecast = PP(5,:)';  % 50% percentile - median prediction
forecast_cases_q75 = PP(6,:)'; % 75% percentile
forecast_cases_q90 = PP(7,:)'; % 90% percentile 
forecast_cases_q95 = PP(8,:)'; % 95% percentile 
forecast_cases_q97p5 = PP(9,:)'; % 97.5% percentile 



indf_ini=665+52+52; % time index of the EW 41 2024
indf_end=716+52+52; % time index of the EW 40 2025

epiweek=[(202441:202452)';(202501:202540)'];  % epidemic weeks to be forecast 


median_cases=mean_cases_forecast(1+(1:52)); % median forecast from EW 41 2024 to EW 40 2025


q2p5=forecast_cases_q2p5(1+(1:52)); % 2.5% percentile  forecast
q5=forecast_cases_q5(1+(1:52)); % 5% percentile  forecast
q10=forecast_cases_q10(1+(1:52)); % 10% percentile  forecast
q25=forecast_cases_q25(1+(1:52));  % 25% percentile  forecast 
q75=forecast_cases_q75(1+(1:52));  % 75% percentile  forecast 
q90=forecast_cases_q90(1+(1:52));     % 90% percentile  forecast
q95=forecast_cases_q95(1+(1:52));     % 95% percentile forecast  
q97p5=forecast_cases_q97p5(1+(1:52));     % 97.5% percentile  forecast


LB95=q2p5;
UB95=q97p5;
LB90=q5;
UB90=q95;
LB80=q10;
UB80=q90;
LB50=q25;
UB50=q75;

T = table(epiweek,median_cases,LB95,UB95,LB90,UB90,LB80,UB80,LB50,UB50);



writetable(T,['..\planilhas\validation_3_Surge_Model_',UF,'.csv'],'Delimiter',',')


range2=indf_ini:min(indf_end,length(dcases_orig));
cases=dcases_orig(range2);

% Generates plots and save as PDF files
max75=max(q75);
max90=max(q90);
max95=max(q95);
max97p5=max(q97p5);

EW_index=indf_ini:indf_end;

figure
subplot(221)
plot(range2,cases,'linewidth',2); % plot known sequence of cases from EW 41 2024 up to EW 52 of 2024
hold on;
plot(EW_index,median_cases,'r','linewidth',2)
plot(EW_index,q25,'k','linewidth',2)
plot(EW_index,q75,'g','linewidth',2)
legend('observed','median forecast','LB50','UB50')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 50% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max75])

subplot(222)
plot(range2,cases,'linewidth',2); % plot known sequence of cases from EW 41 2024 up to EW 52 of 2024
hold on;
plot(EW_index,median_cases,'r','linewidth',2)
plot(EW_index,q10,'k','linewidth',2)
plot(EW_index,q90,'g','linewidth',2)
legend('observed','median forecast','LB80','UB80')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 80% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max90])

subplot(223)
plot(range2,cases,'linewidth',2); % plot known sequence of cases from EW 41 2024 up to EW 52 of 2024
hold on;
plot(EW_index,median_cases,'r','linewidth',2)
plot(EW_index,q5,'k','linewidth',2)
plot(EW_index,q95,'g','linewidth',2)
legend('observed','median forecast','LB90','UB90')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 90% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max95])

subplot(224)
plot(range2,cases,'linewidth',2); % plot known sequence of cases from EW 41 2024 up to EW 52 of 2024
hold on;
plot(EW_index,median_cases,'r','linewidth',2)
plot(EW_index,q2p5,'k','linewidth',2)
plot(EW_index,q97p5,'g','linewidth',2)
legend('observed','median forecast','LB95','UB95')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 95% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max97p5])

print(['..\plots\Validation_3_SModel_',UF],'-dpdf')

close all



