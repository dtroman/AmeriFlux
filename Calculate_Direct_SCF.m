function [ scf ] = Calculate_Direct_SCF( cospectra, freq, bindex )
%Calculate_Direct_SCF takes cospectra and calculates the direct spectral
%correction factor
%   Inputs from load_EP_FullCospectra and make_RH_bins functions
%   Output is vector of spectral correction factors
%   Pascal Polonik, December, 2018

Nfreq = size(cospectra,1);
Nt = size(cospectra,2);
Nvars = size(bindex,2);
% Number of frequency bins. 45 because it is the number of binned 
% frequencies in the EddyPro binned cospectra (at least for 30 min
% averaging)
Nfb = 45; 
scf = NaN(size(bindex));

% Figure out points to bin around (similar to 'binned frequency')
% frequencies in the EddyPro binned cospectra (at least for 30 min
% averaging)
i=int16(logspace(log10(1),log10(Nfreq),Nfb));
idx = @(v) [1 diff(int16(v))];
while any(idx(i)==0)
    fi = find(idx(i)==0);
    i(fi(1):(fi(1)+floor(Nfb/2)))=i(fi(1):(fi(1)+floor(Nfb/2)))+1;
end  

for vi = 1:Nvars
    %figure();    
    bindex_cur = bindex(:,vi);
    % (Nbins below) Minus one because 0 is for points that do not fall in a bin
    Nbins = length(unique(bindex_cur))-1; 
    for b = 1:Nbins
        bi = bindex_cur==b;
        gasspec = nanmedian(cospectra(:,bi,vi),2);      % Gas cospectrum
        Hspec   = nanmedian(cospectra(:,bi,Nvars+1),2); % H cospectrum

        fbinmin = [0, i(1:(end-1))+diff(i)/2];
        fbinmax = [i(1:(end-1))+diff(i)/2, Nfreq+1];
        allidx = 1:Nfreq;

        % Frequency binning (averaging) of the two cospectra
        gasspec_bin = NaN(1,Nfb);
        Hspec_bin = NaN(1,Nfb);
        for fb = 1:length(fbinmin)
            fbi = (allidx>=fbinmin(fb)) & (allidx<fbinmax(fb));
            
            freq_bin(fb) = median(freq(fbi));
            gasspec_bin(fb) = mean(gasspec(fbi));
            Hspec_bin(fb) = mean(Hspec(fbi));
        end
        
        % Calculate transfer function
        tfbin = gasspec_bin./Hspec_bin; % Fratini et al., 2012 Eq. 1
        
        % Interpolate tranfer function back to high frequencies
        tfinterp = interp1(double(i),tfbin,double(1:Nfreq));        
        
        % Limit analysis to high frequencies and remove negative points
        % If there is aliasing in the copectra, this is where changes are required
        tfinterp(tfinterp<0) = min(tfinterp(tfinterp>0));
        tfinterp(1:find(tfinterp>=1,1,'last'))=1;
        
        tfinterp_stack = repmat(tfinterp',1,Nt);

        origspec = cospectra(:,:,Nvars+1);
        
        df = freq(2)-freq(1); % Constant frequency steps
        
        % Equation 3 of Fratini et al., 2012 (no square root)
        Fl_numerator   = sum(origspec.*df./repmat(freq,Nt,1)',1);
        Fl_denominator = sum(origspec.*df./repmat(freq,Nt,1)'.*tfinterp_stack,1);        

        scf_all = Fl_numerator./Fl_denominator;
        scf(bi,vi) = scf_all(bi);
        figure()
        semilogx(tfinterp)
    end
end
% Remove unrealistic corrections in case there are any
scf(scf<0) = NaN;
scf(scf>50)= NaN;
end


