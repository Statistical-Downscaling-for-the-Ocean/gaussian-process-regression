% Script adapted from "get_params_fitgrp_lineP.m" as a test/sample for the
% Line P hackathon to happen on Dec.2025.
% 
% The script reads preprocessed timeseries files (produced in python, test
% is DIC_UMOL_KG, P26, 10 m) obtained from vertically interpolated profiles 
% and fits seasonal cycle using fitgrp and taking into account a linear trend.  
% 
% The test .nc timeseries file provided was produced with: "qc_observed_timeseries.ipybn"
%
% Input needed: 
% Station
% Variable
% 
% History
% 7.Oct.2025: created by Ana C. Franco (ana.franco@bsc.es), based on
% get_params_fitgrp_lineP.m

%% === Initialize ===
clear all
close all
clc

% === User defined ===

% MAKE SURE YOU ARE ON THE RIGHT DIRECTORY!!!
major_st = 26;            % Select the station to analyze
var2plot = 'DIC_UMOL_KG';   % As of 26.jan.2024, should work for 
% sDIC33_UMOL_KG and DIC_UMOL_KG, OXYGEN_UMOL_KG, TEMP_COMB_DEG_C, 
% SALT_COMB_PSS78, SIGMA0_KG_M3, NITRATE_NITRITE_UMOL_KG,
% PHOSPHATE_UMOL_KG, SILICATE_UMOL_KG, PCO2_UATM, ARAG_SAT, PH_TOT, TAslon

plot_timeseries = 'true';

%% ==== Load the files ====
% Loading one file for the moment but leaving the loop to read several. In
% production, the timeseries of each individual depth is in an individual 
% .nc file. A profile for each station could be produced. Code can be 
% extended for many stations at once. 

% Specify the depths to read, file for that interpolated depth should exist
layers = [10]; %layers = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200];

%f1=figure(1);

for l=1:size(layers,2)
    ncfile = ['LineP_',var2plot,'_timeseries_1990_2019_', num2str(layers(l)),'m.nc'];

    %% === Treat the data ===
    % Read variable and date, convert date to matlab date

    ydata_tmp = ncread(ncfile, 'INTERP_VAR_mean'); % DIC timeseries data
    year_tmp = ncread(ncfile, 'YEAR_UTC');         
    month_tmp = ncread(ncfile, 'MONTH_UTC');
    day_tmp = ncread(ncfile, 'DAY_UTC');
    date_tmp = datetime(year_tmp,month_tmp,day_tmp);% Construct the date timeseries

    % Remove duplicates that are actually the mean of two values (comes from the pre-processing). Find them with duplicate index.  
    ydata_idx = ncread(ncfile, 'duplicate_index');
    unique_index = unique(ydata_idx);

    % Keep only the first occurrence of each repeated value
    ydata = [];
    tdata = [];
    for i = 1:size(unique_index,1)
        firstIndex = find(ydata_idx == unique_index(i), 1, 'first');
        ydata = [ydata, ydata_tmp(firstIndex)];
        tdata = [tdata, date_tmp(firstIndex)];
    end

    % Convert date to matlab date and transpose the matrix so that the gpm will
    % recognise it.
    tdata = datenum(tdata)';
    ydata = ydata';

    %% === Prepare input for fitgrp ===
    % 1. Define function handle
    h1 = @(u) [ones(size(u)), u , cos(2*pi*u/(365.25)), sin(2*pi*u/(365.25))];

    % 2. Provide initial guess for parameters. Generate them
    % in relation to the data in this order: constant, slope, amplitude of the
    % cosine and amplitude of the sine terms.
    beta01 = [mean(ydata) (max(ydata)-min(ydata))/(max(tdata)-min(tdata)) std(ydata) std(ydata)];

    % 3. Calculate the gpm. The 'kernel parameters' are 300 = correlation length
    % or link parameter.
    % 0.001 = amplitude of the kernel paremeter, has to be many orders of
    % magnitude smaller than Y, so that it does not "inform" the posterior
    % distribution.
    gpm1 = fitrgp(tdata,ydata,'basisfunction',h1,'beta',beta01,'kernelfunction','squaredexponential','KernelParameters',[300 0.001]);

    % 4. Read beta and sigma
    beta1 = gpm1.Beta; % This has the size 4,1, and values: 
    % 1228.125 0.00107 6.373 31.1206 (for this example 7.Oct.2025, DIC, P26, 10m).
    % Estimated coefficients for the explicit basis functions, stored as a vector. You can define the explicit basis function by using the BasisFunction name-value pair argument in fitrgp.
    
    sy1 = gpm1.Sigma;  % This has just one element with value  11.6860633090684 (for this example).
    % Sigma is the estimated standard deviation of the observation noise, obtained from the GPR model, and stored
    % as a single scalar value: https://es.mathworks.com/help/stats/regressiongp.html

    % H is the set of basis functions evaluated at all training points (i.e., at all available times)
    H1 = h1(tdata)';

    % 5. Calculate mean and covariance matrices
    % This has to do with how to obtain the CI. Betabar1 and beta1 have to be
    % the same. According to Adam, this proves that the theory (the formula in
    % betabar1) agrees with what matlab calculates.
    % betacov1 is the matrix of covariances. The diagonal is the variance of
    % each parameter, and it is used to calculate the credible interval and
    % generate random sets of paramters based on the mean and covariance
    % matrix.
    betabar1 = inv(H1*H1')*H1*ydata; % This should match beta1 and it does, the theory matches... it has a 4x1 shape
    
    betacov1 = (sy1^2)*inv(H1*H1');  % Size 4 by 4 
    % This gives a warning that I have not dealt with yet: Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  3.137945e-17. 
   
    %% === Calculate statistics for each parameter ===

    % 1. Calculate the means of the parameters and store them in a
    % structure called gpr.
    gpr.intercept_mean(l,1) = beta1(1)+beta1(2)*tdata(1);      % Mean intercept
    gpr.slope_yr_mean(l,1) = beta1(2)*365.25;                  % Mean annual slope
    gpr.amplitude_mean(l,1) = sqrt(beta1(3)^2 + beta1(4)^2);   % Mean amplitude
    
    % The variance of each parameter is the diagonal of the betacov1 matrix.
    % 2. Get the standard deviations
    % Intercept
    gpr.intercept_std(l,1) = sqrt(betacov1(1,1));              % Std dev intercept
    
    gpr.slope_std(l,1) = sqrt(betacov1(2,2))*365.25;           % Std dev annual slope

    % The Std dev amplitude was calculated differently. 7.Oct.2025: Could elaborate on this if necesary.   
    % We debated for some time on how to calculate the standard deviation of the amplitude?

    %% === Save the original tdata, ydata, beta1 and betacov1 for every layer to use in python ===
    %layerName = sprintf('gpr.layer%dm', layers(l)); % Assuming layers are named like layer10m, layer20m, etc.
    %eval([layerName '.beta1 = beta1;']);
    %eval([layerName '.betacov1 = betacov1;']);
    %eval([layerName '.tdata = tdata;']);
    %eval([layerName '.ydata = ydata;']);
    %gpr.layers=layers;

    % save params to use later. 
    %fname = ['gpr_params_', var2plot,'_P', num2str(major_st),'_0_200m.mat'];

    % ======= Some other operations/test plots made in the loop =======
    
    % 6. Evaluate function with the mean parameters (beta1) to generate the mean
    % response and slope
    tfit = [tdata(1):1:tdata(end)]';    % The function will be evaluated at all these dates (including gaps). Starts from the first observed date
    yfit = beta1(1) + beta1(2)*tfit + beta1(3)*cos(2*pi*tfit/(365.25)) + beta1(4)*sin(2*pi*tfit/(365.25));
    mean_slope = beta1(1) + beta1(2)*tfit;

    % 7. Evaluate function only for the available dates, to match the timing of observed data and calculate residuals
    yfit_sub = beta1(1) + beta1(2)*tdata + beta1(3)*cos(2*pi*tdata/(365.25)) + beta1(4)*sin(2*pi*tdata/(365.25));


    % 3. Generate 1000 sets of (a0,a1,a2,a3) or (a0,a1) based on the
    % multivariate sampling distribution of these parameters.
    n = 1000; % I typically use 1000  
    R = mvnrnd(beta1,betacov1,n);
    beta3 = R(:,3); beta4 = R(:,4);

    % Plot a PDF for each parameter
    figure;
    titles = {'Beta1', 'Beta2', 'Beta3', 'Beta4'}
    for s = 1:4 % do one plot for each of the 4 columns (beta) 
        subplot(2,2,s);
        ksdensity(R(:,s));
        title(titles(s))
    end

    % Generate 1000 sampels of amplitude A
    sampled_amplitudes = sqrt(beta3.^2 + beta4.^2);
    q025 = quantile(sampled_amplitudes, 0.025);
    q975 = quantile(sampled_amplitudes, 0.975);
    figure; ksdensity(sampled_amplitudes)
    hold on
    xline(q025); xline(q975); 
    hold on
    xline(gpr.amplitude_mean(l,1),'--')
    title('1000 Amplitudes')

    if strcmp(plot_timeseries,'true')
        % Now generate the fitted model for each of the 1000 distributions
        figure
        for r = 1:n
            % Whole fit
            yfit_tmp(r,:) = R(r,1) + R(r,2)*tfit + R(r,3)*cos(2*pi*tfit/(365.25)) + R(r,4)*sin(2*pi*tfit/(365.25));
            hold on
            plot(tfit',yfit_tmp(r,:),'color',[0.8,0.8,0.8])

            % Only slope
            slope_tmp(r,:) =  R(r,1) + R(r,2)*tfit;
            hold on
            plot(tfit',slope_tmp(r,:),'color',[0.8,0.8,0.8])
        end

        % plot the data and the fit
        hold on
        plot(tfit,yfit,'r-','linewidth',1);          % Add the mean function
        hold on
        plot(tfit,mean_slope,'r','linewidth',1);    % Add the mean slope
        hold on
        plot(tdata,ydata,'bo')                      % Plot the original data
        datetick('x', 'mmmyy');

       
% 
%       if save_timeseries_plots == true
%           % Need to find directory to save time series
%       end
     end
% 
%     %% === Calculate the residuals ===
% 
%     res_std(l) = std((ydata - yfit_sub),1);
%     seas_std(l) = nanstd(yfit);
%     total_std(l) = nanstd(ydata);
% layers(l)
end

%% === Plot the amplitude vertical profile ===
% 
% f2 = figure(2);
% f2 = fill([gpr.amplitude_mean(:)-gpr.amplitude_std(:);flipud(gpr.amplitude_mean(:)+gpr.amplitude_std(:))],[layers';flipud(layers')],[.7,.7,.7],'linestyle','none','LineWidth',0.2);
% hold on
% h1 = plot(gpr.amplitude_mean(:),layers,'k.-','LineWidth',2); % Plot the amplitude
% set(gca,'YDir','reverse')
% hold on
% h2 = plot(res_std(:),layers,'r.-','LineWidth',2); % Plot the standard deviation of residuals
% hold on
% h3 = plot(total_std(:),layers,'b.-','LineWidth',2); % Plot the total standard deviation
% hold on
% h4 = plot(seas_std(:),layers,'g.-','LineWidth',2);
% %     set(gca,'Ydir','reverse')
% if strcmp(var2plot,'OXYGEN_UMOL_KG')
%     xlim([0 40]); xlabel('O2, µmol/kg')
% elseif strcmp(var2plot,'DIC_UMOL_KG')
%     xlim([0 50]); xlabel('DIC, µmol/kg')
% elseif strcmp(var2plot,'sDIC33_UMOL_KG')
%     xlim([0 40]); xlabel('sDIC33, µmol/kg')
% elseif strcmp(var2plot,'temp')
%     xlim([0 5]); xlabel('Temp, degC')
% elseif strcmp(var2plot,'sigma0')
%     xlim([0 1.2]); xlabel('Sigma0, kgm-3')
% elseif strcmp(var2plot,'pH_2mlr')
%     xlim([0 0.03]); xlabel('pH 2mlr')
% elseif strcmp(var2plot,'arag_2mlr') || strcmp(var2plot,'aragonite')
%     xlim([0 0.4]); xlabel('arag 2mlr')
% elseif strcmp(var2plot,'salt')
%     xlim([0 0.5]); xlabel('salt')
% elseif strcmp(var2plot,'NO3')
%     xlim([0 5]); xlabel('NO3 µmol/kg')
% end
% ylim([10 100])
% set(gca,'FontSize',16)
% ylabel('Depth (m)')
% legend([h1 h2 h3 h4],'Seasonal amplitude','std dev residuals','total std dev', 'seas std')
% grid on
% 

