function [ bindex, time, data, header ] = make_RH_bins( path,start_date,end_date,bin_opt )
%make_RH_bins takes full EddyPro outputs and returns a vector of bin numbers
%   path is the path and file name of the EddyPro fulloutput file
%   start_date and end_date have the format 'YYYY-MM-DD'
%   bin_opt is 'none' or 'quantile' - could be adapted to include more


% Load full output to read relative humidities (RH)
[data,header,~] = xlsread(path);
units = header(3,3:end);
header = header(2,3:end);
fid = fopen(path);
dates = textscan(fid,'%*s%s%*[^\n]','Delimiter',',','HeaderLines',3);
fclose(fid);
yr = char(dates{1});
yr = str2num(yr(:,1:4));

daten = datenum(yr,0,data(:,strcmp(header,'DOY')));

idx = (daten>datenum(start_date)) & (daten<=datenum(end_date)+1);
data=data(idx,:);
time=daten(idx);
RH = data(:,strcmp(header,'RH'));

% Get the indices for each selected option
bindex = zeros(size(RH,1),numel(bin_opt));
for boi = 1:numel(bin_opt)
    % Apply RH bin options to get bin limits
    switch bin_opt{boi}
        case 'quantile'
            RH_bins = [nanmin(RH)-0.01, quantile(RH,3), nanmax(RH)+0.01];
        case 'none'
            RH_bins = [0,100];
    end

    % Label each point with a bin (e.g. 0=exclude, 1=first bin, etc.)    
    for bi = 1:(length(RH_bins)-1)
        thisbin = (RH>=RH_bins(bi)) & (RH<RH_bins(bi+1));
        bindex(thisbin,boi) = bi;    
    end
    
end
end

