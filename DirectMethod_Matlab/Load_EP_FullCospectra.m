function [spec,spec_time,freq] = Load_EP_FullCospectra(path, start_day, end_day, variable)
%LoadDavisSpectra loads the Ameriflux Davis spectral outputs from EddyPro
%   Start Day and End Day are input as date strings ('yyyy-mm-dd')
%   Path is the path from the current folder to get to folder with all the
%       Path is often: '\\arm.lbl.gov\davis\alldaycsv'
%   data (organized by date)
%       Options: 'spec','cospec','both'    
%
%   Pascal Polonik, Dec 2018

% Number of days selected
nday = datenum(end_day)-datenum(start_day)+1;
if nday <= 0
    warning('Warning: End day is before start day')
end

Nvars = numel(variable);

%possible_headers = {'natural_frequency','normalized_frequency','f_nat*cospec(w_ts)','f_nat*cospec(w_co2)','f_nat*cospec(w_h2o)'};

f = dir([path '*.csv']);
fnames = cell(numel(f),1);
for fi = 1:numel(f)
    fnames{fi} = f(fi).name;
end

% Read first file to get info
[~, ~, ~, meta1] = read_cospectrum(path,f(1));
Hz = meta1(1);
avg_period = meta1(4);
nseg = 24*60/avg_period;  % Number of averaging periods in a day
ppf = 2^floor(log2(avg_period*60*Hz/2)); % Points per file

freq = 0:Hz/2/ppf:Hz/2;
freq = freq(2:end);

%spec=[];
spec = NaN(ppf,nday*(24*60/avg_period),Nvars);
spec_time=[];

day = start_day;
day = datestr(datenum(day),29);
tct=0;
for d = 1:nday
    % Run through each segment of the day (e.g. 48 half hours)
    for h = 1:nseg
        % Loop goes from first non-zero time to midnight of the next day
        % Therefore, before the last time, one day needs to be added
        tct=tct+1;
        if h==nseg
            day = datestr(datenum(day)+1,29);
        end
        
        hstr = datestr(datenum(day)+(h*avg_period)/(24*60),'HHMM');

        daystr = strjoin(strsplit(day,'-'),'');
        disp(['Loading... ' day ' ' hstr])

        % See if file exists
        matchi = strmatch([daystr '-' hstr],fnames);
        
        % If it exists, continue loading. Otherwise make NaNs
        if ~isempty(matchi)
            [a,b,c] = fileparts(f(matchi).name);     %a= folder, b= file, c= ext
            t = regexp(a,'\','split');     
            [spec_day, spec_time_day, header_day, ~] = read_cospectrum(path, f(matchi));          

            % xlsread messes up the finding of the data, so remove the first two rows
            if size(spec_day,1) == ppf + 2
                spec_day = spec_day(3:end,:);
            end
            spec_day_use = NaN(ppf,Nvars);
            
            for vi = 1:Nvars
                gasheader = ['f_nat*cospec(w_' variable{vi} ')'];

                % Run through each existing column and fill spec_day_use with data from spec_day
                for col = 1:length(header_day)
                    hi = strcmp(gasheader,header_day{col}); 
                    if sum(hi)==1 % Make sure something reasonable is found (somteimes NaN or empty header_day)
                        try
                            spec_day_use(:,vi) = spec_day(:,col);                            
                        catch
                            error(['Bad file: ' daystr '-' hstr])
                        end
                    end
                end
            end
            
            % Append selected cospectra
            spec(:,tct,:) = spec_day_use;
            spec_time = [spec_time; spec_time_day];
         
        else
            % What to do if no file is found
            spec_time = [spec_time; NaN];
        end
                
    end
    
end

end

