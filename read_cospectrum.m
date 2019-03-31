function [spec,timeseries,header,meta] = read_cospectrum(path,d)

% batch loader for EddyPro binned (co)spectra results
% inputs:
% path = directory where files are located
% d = structure from DIR with all files to load
% Assumes that there is a date in the file name with the format
%   yyyymmdd-HHMM
%
% outputs:
% spec = stacked output of all frequency bins x columns of (co)spectra
% timeseries = datenum of each file used
% header = column names
% meta = [sampling frequency,sample height, wind speed, averaging period]


% Oringinally Stephen Chan, Oct 2014
% Revised Pascal Polonik, Dec 2018

%%

spec = [];
for i = 1:numel(d)    
    filename = [path d(i).name];
    f=fopen(filename);
    header_all = textscan(f, '%s', 13, 'delimiter', '\n');   
    Hz = strsplit(char(header_all{1}{6}),'_');
    Hz = str2double(Hz(end));
    height = strsplit(char(header_all{1}{7}),'_');
    height = str2double(height(end));
    ws = strsplit(char(header_all{1}{8}),'_');
    ws = str2double(ws(end));
    avg_period = strsplit(char(header_all{1}{9}),'_');
    avg_period = str2double(avg_period(end));
    meta = [Hz,height,ws,avg_period];
    fclose(f);
    header = strsplit(char(header_all{1}{13}),',');
    dat = csvread(filename,13,0);
    dstringi = regexp(d(i).name,'\d{8}-\d{4}');
    timeseries(i,1) = datenum(d(i).name(dstringi:(dstringi+12)),'yyyymmdd-HHMM');     
    spec = [spec; dat];
end

% convert -9999 to NaN
spec(spec==-9999) = NaN;
