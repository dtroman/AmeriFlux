% Setup to calculate Direct spectral correction method from EddyPro full
% cospectra as decribed in Polonik et al., 2019 Ag & Forest Met
% Script is publically available. If used, please include citation.
% Created by Pascal Polonik and Stephen Chan, Dec 2018

% Please be aware that all cospectra in the selected time are loaded 
% simultaneously and therefore require substantial RAM and time. 
% Cospectra are saved as .mat files when running the first execution.

%% --------------------------------------------------------------------- %%
% OPTIONS
%%-----------------------------------------------------------------------%%

% Set path to folder containing full EddyPro cospectra
% Make sure it has a slash at the end
path_cospec = '\Path\To\Cospectra\';


% Set path to file containing full half-hour output (contains RH)
path_fullout = '\Path\To\Fullout\file.csv';

% Name of sensor (used for saving file)
sensorName = 'SensorNameSave';

% Select Gases
% Cell of strings. Current options: 'co2','h2o'
% No need to include 'ts' because it is added later.
gases = {'co2','h2o'};

% Select relative humidity bin type
% Current options are 'none' and 'quantile', but initialization below could
% easily be altered to construct different bins
% Cell of strings, same size as gases.
bin_opt = {'none','quantile'};

% Set start and end dates ('YYYY-MM-DD')
start_date = '2015-07-05'; % Inclusive: starts at 00:00 of given day
end_date = '2015-10-05'; % Inclusive: ends at 00:00 of the following day

% File to save out cospectra, because loading each file takes several mins.
% If new date range is selected, cospectra will be re-loaded and saved
% separately 
matout = [sensorName '_DirectCospec_' start_date '_' end_date '.mat'];

%% --------------------------------------------------------------------- %%
% SANITY CHECKS
%%-----------------------------------------------------------------------%%
if ~(numel(gases)==numel(bin_opt))
    error('Must choose one bin option per selected gas!')
end
for i=1:numel(gases)
    if ~(strcmp(gases{i},'co2') || strcmp(gases{i},'h2o') || strcmp(gases{i},'ch4'))
        error('Selected gases must be co2, h2o, or ch4')
    end
    if ~(strcmp(bin_opt{i},'none') || strcmp(bin_opt{i},'quantile'))
        error('Selected bin_opt must be either none or quantile')
    end
end
if datenum(end_date)<datenum(start_date)
    error('End date must be after start date!')
end

%% --------------------------------------------------------------------- %%
% RUN SPECTRAL CALCULATIONS
%%-----------------------------------------------------------------------%%

% Variables to load: gases and ts
vars = gases;
vars{numel(vars)+1}='ts';

% Load cospectra. This can be slow depending on the time range.
% Often requires quite a bit of RAM since all cospectra are loaded simultaneously
% Returns co2, h2o, ts cospectra and the corresponding frequencies
% Each cospec variable has a size of Nfreq x Ntimes with each bin in a cell

if exist(matout,'file')
    load(matout)
    disp(['Loaded cospectra from ' matout])
else
    [cospec,time,freq] = Load_EP_FullCospectra(path_cospec, start_date, end_date, vars);
    save(matout,'cospec','time','freq')
    disp(['Saved cospectra to ' matout])
end

% Insert index here to apply filters
% Set points to NaN that should not be used (e.g. cospec(:,index,:) = NaN)
% Calculate Direct spectral correction

[bindex,btime,data,header] = make_RH_bins(path_fullout,start_date,end_date,bin_opt);

cospec_indexed = cospec;
%Example when index is used:
%index = abs(data(:,(strcmp(header,'wind_dir')))-180)<60;
%cospec_indexed(:,~index,:)=NaN;

direct_scf = Calculate_Direct_SCF(cospec_indexed,freq,bindex);

%% --------------------------------------------------------------------- %%
% PLOT
%%-----------------------------------------------------------------------%%

figure(1);
clf
tmp = cospec_indexed;

for i = 1:numel(vars)
    semilogx(freq, nanmedian(tmp(:,:,i),2));
    hold on
end
set(gca,'XScale','Log')
legend(vars,'FontSize',14)

figure(2);
clf
plot(time,direct_scf,'.')
datetick('x')
ylabel('Direct Correction','FontSize',16)
legend(gases)
ylim([-0.5,2])


