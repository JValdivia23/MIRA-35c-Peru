function fawa = get_wetatt(npw1,time,met)
% get_wetatt - Wet antenna attenuation retrieval
%     
%     This function returns a factor of wet antenna attenuation.
%     The basis is on radiometric noise power, which increase in the rain
%     fall.
%       fawa = get_wetatt(npw1,time,method)
%     Where, fawa is the factor of wet antenna attenuation, and, npw1, is
%     the radiometric noise. We ca chose linear o exponential method (1 or 2)
% 
% Created by: Jairo Valdivia                                    Jan - 2017

% For corrections to wet antenna attenuation: (from WetAtt)
% % mm_p(end) and mm_r(end) from testingDSD.m
% ratt = 1.439; %cum of rain att
% rref = 7.844; % cum reference (rain gauge)
% F = rref/ratt; % factor of wet att
% 
% npw1=[chunk1.npw1,chunk2.npw1];
% npwt=npw1(1:2320); % just until rain fall
% ref = 85;
% nwr1 = sum(npwt(npwt>ref+1))/numel(npwt(npwt>ref+1)); %mean of noise with rain, +1 to noise ref
% mwnRef = nwr1 - ref; % mean wet noise with noise reference
% 
% cwa = F/mwnRef  % constant of wet attenuation --> DSD.m

t = time-round(time);%(0:1439)*1/1440;
smth=round(60/86400/mean(diff(time)));
a0=87.2594; a1=-3.0484; b1=-1.4592; a2=-1.8176; b2= 1.0342;% a3=-0.0045; b3=0.1425; %92.5
% a4=0.1996; b4=0.7347;
% a0=90; a1=2.8772; b1=-1.6040; a2=-0.0257; b2=1.5017; a3=0.2606; b3=0.6029;
% a4=-0.0833; b4=0.6677;
X=a0+a1.*sin(2*pi*t)+b1.*cos(2*pi*t)+...
    a2.*sin(2*pi*2*t)+b2.*cos(2*pi*2*t);%+...
%     a3.*sin(2*pi*3*t)+b3.*cos(2*pi*3*t);%+...
%     a4.*sin(2*pi*4*t)+b4.*cos(2*pi*4*t);

if ~exist('met','var'), met=1; end
switch met
    case 1
        
    cwa = 1; %0.5 Constant of wet att from (WetAtt.m)
     % needed 86 to make better
    smth=10; % smoothing data
    fawa = NaN(size(npw1));
    for nd = 1:numel(npw1)
        ref = X(nd);
        if nd == 1, ii = 1; sm_fawa = []; end
       if ii > smth, ii=1; end
       if npw1(nd)>ref+1/cwa
            sm_fawa(ii) = cwa*(npw1(nd)-ref);                
       else
           sm_fawa(ii) = 1;
       end
       fawa(nd) = mean(sm_fawa);
        ii=ii+1;
    end
    
 % para proporción exponencial: A = exp(kdN)
    case 2
        
    cwa = 1/20; % .73  15.0626 25.1817 
%     smth=10; % smooth of data
    fawa = NaN(size(npw1));
    for nd = 1:numel(npw1)
       ref = X(nd);
        if nd == 1, ii = 1; sm_fawa = []; end
       if ii > smth, ii=1; end
       if npw1(nd)>ref
            sm_fawa(ii) = 10^(cwa*(npw1(nd)-ref));                
       else
           sm_fawa(ii) = 1;
       end
       fawa(nd) = mean(sm_fawa);
        ii=ii+1;
    end
    % para proporción exponencial: A = dN^{1/k}
    case 3
        
    cwa =   1.34;%1.8;%0.78; % 
    kq  =   1.1;%0.94;%
%     smth=10; % smooth of data
    fawa = NaN(size(npw1));
    for nd = 1:numel(npw1)
       ref = X(nd);%89;
        if nd == 1, ii = 1; sm_fawa = []; end
       if ii > smth, ii=1; end
       if npw1(nd)>ref+(1/cwa)^kq
            sm_fawa(ii) = cwa*(npw1(nd)-ref)^(1/kq);                
       else
           sm_fawa(ii) = 1;
       end
       fawa(nd) = mean(sm_fawa);
        ii=ii+1;
    end
end

% 
%     
%         for nd = 1:numel(npw1)
%         if nd == 1, ii = 1; sm_fawa = []; end
%         if nd <= smth
%             if npw1(nd)>ref
%                 fawa(nd) = 10^(cwa*(npw1(nd)-ref));
%                 sm_fawa = [sm_fawa, fawa(nd)];
%                 fawa(nd) = mean(sm_fawa);
%             else
%                 sm_fawa = [sm_fawa, 1];
%                 fawa(nd) = mean(sm_fawa);            
%             end                   
%         end
% 
%         if nd > smth
%            if ii > smth, ii=1; end
%            if npw1(nd)>ref
%                 sm_fawa(ii) = 10^(cwa*(npw1(nd)-ref));
%                 fawa(nd) = mean(sm_fawa);
%            else
%                fawa(nd) = mean(sm_fawa);
%            end
%             ii=ii+1;
%         end
%     end