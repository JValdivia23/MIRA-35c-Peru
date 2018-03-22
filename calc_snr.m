function [snr, vlx, vvx]=calc_snr(spc,noise)
% Using the signal (spc) and noise (noise), this function
% calculates frist and second moment of Doppler Spectra (vlx & vvx)
% also calculates signal to noise (snr)
% By:	Jairo Valdivia

xx = find(spc>0);
sig = spc(xx);%-hild;
pwx = sum(sig);
if ~isempty(xx)
    vlx = sum(spc(xx).*xx)/pwx;
    vvx = sum(spc(xx).*xx.*xx)/pwx;
    snr=NaN(size(spc));
    if length(xx)> 2
        snr(xx)=sig/noise;% vlx=vlx; 
        vvx=sqrt(vvx - vlx.*vlx);
    return    
    else
        snr(xx)=sig/noise;% vlx=vlx; 
        vvx=NaN;
    end
else
    snr=NaN(size(spc)); vlx=NaN; vvx=NaN;
    return
end