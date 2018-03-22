function mietable = mietab
% This script generates a structure with Mie factor corrections and single
% particle backscattering cross sections for diferentes diameters and temp
% tic
if exist('mietable.mat','file')
    disp('Mietable already exist')
    return ; 
end

xmt = 34.85*10^9; %frecuency of MIRA 35c 
c_light = 2.99792458*10^8;
lambda = c_light/xmt;
%  tabulando n y k con w1=0.62cm y w2=1.24cm
%   this values has been exacted from Oklahoma University web page:
%   http://www.ou.edu/radar/
% w1=0.0062; w2=0.0124;
% w1n=[3.10,3.45,3.94,4.44]; w1k=[1.77,2.04,2.37,2.59];
% w2n=[4.15,4.75,5.45,6.15]; w2k=[2.55,2.77,2.9,2.86];
% w1kwq=[0.7921,0.8312,0.8726,0.8926];
% w2kwq=[0.8902,0.9055,0.9152,0.9193];
% 
% brit_n_n = ((lambda-w1).*w2n + (w2-lambda).*w1n)/(w2-w1);
% brit_n_k = ((lambda-w1).*w2k + (w2-lambda).*w1k)/(w2-w1);
% brit_temp = [-8,0,10,20];
% brit_kwq = ((lambda-w1).*w2kwq + (w2-lambda).*w1kwq)/(w2-w1);% constante |k|^2
complexindex % Complex index generator
bri_temp = temp1; bri_nt = numel(bri_temp);%linspace(-20,35,56); bri_nt=length(bri_temp); 
bri_n_n = val_n;%interp_jairo(bri_temp,brit_temp,brit_n_n);
bri_n_k = val_k;%interp_jairo(bri_temp,brit_temp,brit_n_k);
bri_d = 0.00001:0.00001:0.008; bri_nD = numel(bri_d); %drop diamter in m
bri_kwq = val_kwq;%interp_jairo(bri_temp,brit_temp,brit_kwq);


% mietable = NaN(bri_nD,bri_nt);
sbcross = NaN(bri_nD,bri_nt);
raycross = NaN(bri_nD,bri_nt);
extinc = NaN(bri_nD,bri_nt);
mietable = struct('Drop_diameter',bri_d*1000,'Temperture',bri_temp,...
    'SBcross',sbcross,'RayCross',raycross,'Extinction',extinc);

kwq = 0.93;     % constante |k|^2 del agua utilizada para calcular la constate de radar

for it = 1:bri_nt
    nref = complex(bri_n_n(it),bri_n_k(it));
    for id = 1:bri_nD
        [s1,s2,qext,qsca,qback,gsca] = bhmie(pi*bri_d(id)/lambda,nref,1);
        %mietab.mietable(id,it) = qback*lambda^4/(kwq*bri_d(id)^2*pi^3);
        mietable.SBcross(id,it) = qback*pi^2*(bri_d(id))^2; %pi/4*D;
        mietable.RayCross(id,it) = bri_kwq(it)*bri_d(id)^6*pi^5/lambda^4;
        mietable.Extinction(id,it) = qext*pi*(bri_d(id)/2)^2;
    end
end
save('mietable.mat','mietable')
% toc
end

