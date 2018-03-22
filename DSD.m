% TO CALCULATE DROP SIZE DRISTRIBUTION
function [dropsize, chunk] = DSD(filename)
% DSD - Drop Size Distribution retrieval
%
%   This function use the Doppler specta to retrieve the DSD and others
%   microphisical paremeters of rainfall. 
%       [dropsize, chunk] = DSD(filename)
%   where 'filename' are zspc files and 'dropsize' and 'chunk' are
%   structures with parameters and spectra.
%
%   UPDATES:
%   01/19/2018 - added PIA algorithm/DSD now is drop density (b. drop concentration)

if ~exist('filename','var'), filename = '20151229_0000.zspc'; end
try 
gunzip(filename)
    p=find(filename=='.'); p=p(end);
%     filename2=filename;
    filename=filename(1:p-1);
    chunk = spc_read2(filename);
    delete(filename);
catch
    error('Invalid data or no compresed')
end


disp('Getting microphysical parameters...')

if exist('mietable.mat','file'), load('mietable.mat') ;
else
    mietable=mietab;
end


f_co = chunk.Co_Spc_Mtr;
UTC = chunk.UTC;
COFA = chunk.COFA_co;
RC = chunk.RadarConst5;
npw1 = chunk.npw1;
range = chunk.range;
dwells=length(UTC);
nbins = chunk.Process_Param.sft;
nrange = length(range);
kwq=0.89; %0.93 indice complejo de refracion estandar | a 10°c y 34.85 Ghz kwq = 0.89
c = 299792458;
xmt = 34.85*10^9; %frecuencia
lambda = c/xmt; %wavelength

% For nyquist correction:
fix=20; %~3 m/s
% shift = nbins/2+fix; %bin shift for corr 
shift = fix-1; % (IN NEW FILES: shift=fix-1 )
f_co=fliplr(f_co);
f_co=[f_co(:,nbins-shift+1:end,:) f_co(:,1:nbins-shift,:)];
  % 
HSDV_co=chunk.HSDV_co;
ny_vel = c * chunk.Process_Param.prf / (4.0*xmt);
s_vel=fix*ny_vel*2/nbins;
vel = 2*ny_vel*((1:nbins)-nbins/2)/nbins;
vel = vel - s_vel;  %vel(-vel<0)=NaN;
masl=3230; % meters above ground level
nue = 1+3.68*10^-5.*(range+masl)+1.71*10^-9.*(range+masl).^2; 

% fawa = get_wetatt(npw1,UTC,3); % factor of wet antenna attenuation get_wetatt(noise, met)
%% Eta

eta = NaN(nrange,nbins,dwells);
vlx = NaN(nrange,dwells); vvx = NaN(nrange,dwells);
for nd = 1:dwells
%     if npw1(nd)>ref+1/cwa, fawa = cwa*(npw1(nd)-ref); else fawa = 1; end %factor of att in wet antenna
    %noise = db2pow(npw1(nd))*noinor1*noinor2;%20.78*npw1(nd)+1831;%HSDV_co(rg,nd)*nbins
    for rg = 1:nrange
        [snr, v, vv]=calc_snr(f_co(rg,:,nd)-HSDV_co(rg,nd),npw1(nd)); % noise/nbins
        vlx(rg,nd)=interp_jairo(v,1:nbins,vel);
        vvx(rg,nd)=2*ny_vel/nbins*vv;
        eta(rg,:,nd) = snr*RC(nd)*COFA(rg,nd)*((range(rg))/5000)^2*10^-18*pi^5*kwq/lambda^4; % fawa(nd)*       
    end
end

Ze=squeeze( sum(eta,2,'omitnan')/(10^-18*pi^5*kwq/lambda^4) );
noval=Ze==0;
Ze(noval)=NaN;
%%
D=NaN(nrange,nbins);
for n = 1:nrange
    velnorm=-vel./nue(n);
    for b = 1:nbins
        if velnorm(b)> 9.6, velnorm(b)=NaN; end %est. is 9.25!
            D(n,b) = ( -1.667 * log( 0.9369 - 0.097087 * (velnorm(b)) ));
            %D(n,b) =  1/0.3*log(10.3./(9.65 -velnorm(b))); %9.65 -> 10.65
        if D(n,b)< 0.1, D(n,b)=NaN; end 
    end
end

% eta_v = eta./(ny_vel/nbins*2);%repmat(vel,[nrange 1 dwells]);
dv_dd = 6.18*exp(-0.6.*D).*repmat(nue,[1 nbins]);
% eta_d = eta_v.*repmat(dv_dd,[1 1 dwells]);
dD = (ny_vel/nbins*2)./dv_dd;

t_sup=10; % 0.6ºC por cada 100 m
t_min=t_sup-6/1000*range(end);
temp=linspace(t_sup,t_min,nrange)';temp=round(temp);
drop=round(mietable.Drop_diameter,2); Drd=round(D,2); %-dD/2
sigma=NaN(nrange,nbins);
sige=NaN(nrange,nbins);

for rg=1:nrange
    t=find(temp(rg)==mietable.Temperture);
    if ~isempty(t),t=t(1);
        for n =1:nbins
            d=find(Drd(rg,n)==drop);
            if ~isempty(d)
                d=d(1);
                sigma(rg,n)=mietable.SBcross(d,t);
                sige(rg,n)=mietable.Extinction(d,t);
            end
        end
    end
end

dsd = eta./repmat(sigma,[1 1 dwells]);
%dsd = eta_d./repmat(sigma,[1 1 dwells]);
wetfact=[];
% Model of error E(R)
RR = squeeze(sum(0.0006.*pi.*repmat(D,[1 1 dwells]).^3.*repmat(-vel,[nrange 1 dwells]).*dsd,2,'omitnan')); RR(RR==0)=NaN;
a = 0.2273 ; b = 1.65; % from fit exp1
mu = 19.2; % Mean HS2
sigma = 0.1777; % DS HS2
ep = 15.8; ga = 45.3;
nref=(10*log10(HSDV_co(end,:)*nbins)-mu); nref(nref<0)=0;
Er = a*exp(b/sigma*(nref)).*(1-exp(-ep.*RR(4,:))).*(1-exp(-ga.*abs(nref)));
% Er(Er<0.01)=0;
R2 = RR(4,:)+Er;
wetfact = R2./RR(4,:); wetfact(isnan(wetfact))=1;
W_ATT = repmat(wetfact,[nrange, 1, nbins]);
W_ATT = permute(W_ATT,[1 3 2]);
dsd = dsd.*W_ATT;
% end of model

RR = squeeze(sum(0.0006.*pi.*repmat(D,[1 1 dwells]).^3.*repmat(-vel,[nrange 1 dwells]).*dsd,2,'omitnan')); RR(RR==0)=NaN;
LWC = squeeze(sum(pi/6*repmat(D,[1 1 dwells]).^3.*dsd,2,'omitnan')); LWC(LWC==0)=NaN;
Z = squeeze(nansum(repmat(D,[1 1 dwells]).^6.*dsd,2));
noval=Z==0;
Z(noval)=NaN;
RR(noval)=NaN;
LWC(noval)=NaN;
[dsd, pia, K]= get_pia(dsd,sige,diff(range(1:2)),4);

dsd=log10(dsd./repmat(dsd,[1 1 length(dwells)]));

% 'dD',dD,
dropsize = struct('time',UTC,'range',range,'vel',vel,'drops',D,'dsd',dsd,...
    'RR',RR,'LWC',LWC,'Z',Z,'dD',dD,...
    'Ze',Ze,'vlx',vlx,'vvx',vvx,'wetatt',wetfact,'fix',fix,'Er',Er,'PIA',pia); %

disp('Complete.')
return

