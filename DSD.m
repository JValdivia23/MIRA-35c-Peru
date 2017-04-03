% TO CALCULATE DROP SIZE DRISTRIBUTION
function [dropsize, chunk] = DSD(filename)
if ~exist('dpath','var'), dpath = 'J:\Otros\Radar\'; end
if ~exist('gpath','var'), gpath = 'J:\Otros\Radar\Plots\Spectra\Co_spc\'; end
if ~exist('filename','var'), filename = '20151229_0000.zspc'; end

% t=ncread([dpath,file,'.mmclx'],'time');
% ndwells=length(t);

chunk = spc_read2(filename);

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
shift = nbins/2+fix; %bin shift for corr 
% shift = fix-1; % (IN NEW FILES: shift=fix-1 )
f_co=fliplr(f_co);
f_co=[f_co(:,nbins-shift+1:end,:) f_co(:,1:nbins-shift,:)];
  % 
HSDV_co=chunk.HSDV_co;
ny_vel = c * chunk.Process_Param.prf / (4.0*xmt);
s_vel=fix*ny_vel*2/nbins;
vel = 2*ny_vel*((1:nbins)-nbins/2)/nbins;
vel = vel - s_vel;  %vel(-vel<0)=NaN;
%estoy quitando el 3300 (0)
masl=3230; %0
nue = 1+3.68*10^-5.*(range+masl)+1.71*10^-9.*(range+masl).^2; 

fawa = get_wetatt(npw1,2); % factor of wet antenna attenuation get_wetatt(noise, met)
%% Eta

eta = NaN(nrange,nbins,dwells);
vlx = NaN(nrange,dwells); vvx = NaN(nrange,dwells);
for nd = 1:dwells
%     if npw1(nd)>ref+1/cwa, fawa = cwa*(npw1(nd)-ref); else fawa = 1; end %factor of att in wet antenna
    %noise = db2pow(npw1(nd))*noinor1*noinor2;%20.78*npw1(nd)+1831;%HSDV_co(rg,nd)*nbins
    for rg = 1:nrange
        [snr, v, vv]=calc_snr(f_co(rg,:,nd)-HSDV_co(rg,nd),npw1(nd)); % noise/nbins
        vlx(rg,nd)=interp_jairo(v,1:nbins,vel);
        vvx(rg,nd)=interp_jairo(vv,1:nbins,vel);
        eta(rg,:,nd) = fawa(nd)*snr*RC(nd)*COFA(rg,nd)*((range(rg))/5000)^2*10^-18*pi^5*kwq/lambda^4; %        
    end
end

Ze=squeeze( sum(eta,2,'omitnan')/(10^-18*pi^5*kwq/lambda^4) );
Ze(Ze==0)=NaN;
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
% dv_dd = 6.18*exp(-0.6.*D).*repmat(nue,[1 nbins]);
% eta_d = eta_v.*repmat(dv_dd,[1 1 dwells]);
% dD = (ny_vel/nbins*2)./dv_dd;

t_sup=10; % 0.6ºC por cada 100 m
t_min=t_sup-6/1000*range(end);
temp=linspace(t_sup,t_min,nrange)';temp=round(temp);
drop=round(mietable.Drop_diameter,2); Drd=round(D,2); %-dD/2
sigma=NaN(nrange,nbins);
for rg=1:nrange
    t=find(temp(rg)==mietable.Temperture);
    if ~isempty(t),t=t(1);
        for n =1:nbins
            d=find(Drd(rg,n)==drop);
            if ~isempty(d)
                d=d(1);
                sigma(rg,n)=mietable.SBcross(d,t);
            end
        end
    end
end

dsd = eta./repmat(sigma,[1 1 dwells]);
%dsd = eta_d./repmat(sigma,[1 1 dwells]);

RR = squeeze(sum(0.0006.*pi.*repmat(D,[1 1 dwells]).^3.*repmat(-vel,[nrange 1 dwells]).*dsd,2,'omitnan'));
LWC = squeeze(sum(pi/6*repmat(D,[1 1 dwells]).^3.*dsd,2,'omitnan'));

% % Z = NaN(nrange,dwells);
% RR = NaN(nrange,dwells);
% % LWC = NaN(nrange,dwells);
% 
% for nd = 1:dwells
%     for rg = 1:nrange
%         % without dD! 'cos dsd [m^-3]  .*dD(rg,:)
%         %Z(rg,nd)=sum(dropsize1.drops(rg,:).^6.*dsd(rg,:,nd).*dD(rg,:),'omitnan');
%         RR(rg,nd)=sum(0.0006.*pi.*dropsize1.drops(rg,:).^3.*-vel.*dsd(rg,:,nd),'omitnan');%5.5
%         %LWC(rg,nd)=sum(pi/6*dropsize1.drops(rg,:).^3.*dsd(rg,:,nd).*dD(rg,:),'omitnan');
%     end
% end

% 'dD',dD,
dropsize = struct('time',UTC,'range',range,'vel',vel,'drops',D,'dsd',dsd,...
    'RR',RR,'LWC',LWC,...
    'Ze',Ze,'vlx',vlx,'vvx',vvx,'fix',fix);
return

