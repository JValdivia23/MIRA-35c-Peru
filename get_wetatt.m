function fawa = get_wetatt(npw1,met)
% get_wetatt - Wet antenna attenuation retrieval
%     
%     This function returns a factor of wet antenna attenuation.
%     The basis is on radiometric noise power, which increase in the rain
%     fall.
%       fawa = get_wetatt(npw1,method)
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
if ~exist('met','var'), met=1; end
switch met
    case 1
        
    cwa = 0.4629; % Constant of wet att from (WetAtt.m)
    ref = 86; % needed 86 to make better
    smth=10; % smooth of data
    fawa = NaN(size(npw1));
    for nd = 1:numel(npw1)
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
        
    cwa = 1/20; %15.0626 25.1817 
    ref = 85; % 
    smth=10; % smooth of data
    fawa = NaN(size(npw1));
    for nd = 1:numel(npw1)
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