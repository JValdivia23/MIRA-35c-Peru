% TO CALCULATE DROP SIZE DRISTRIBUTION
function dropsize = DSD(dpath,filename,ndwells)
if ~exist('dpath','var'), dpath = 'J:\Otros\Radar\'; end
if ~exist('gpath','var'), gpath = 'J:\Otros\Radar\Plots\Spectra\Co_spc\'; end
if ~exist('filename','var'), filename = '20151229_0000.zspc'; end

% t=ncread([dpath,file,'.mmclx'],'time');
% ndwells=length(t);

chunk = spc_read(dpath,filename,ndwells);
mietable=mietab;

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
shift=63-fix; %bin shift for corr
f_co=[f_co(:,nbins-shift+1:end,:) f_co(:,1:nbins-shift,:)];
f_co=fliplr(f_co);
% HSDV_co=chunk.HSDV_co;
ny_vel = c * chunk.Process_Param.prf / (4.0*xmt);
s_vel=fix*ny_vel*2/nbins;
vel = 2*ny_vel*((1:nbins)-nbins/2)/nbins;
vel = vel - s_vel;  %vel(-vel<0)=NaN;
%estoy quitando el 3300 (0)
masl=3230; %0
nue = 1+3.68*10^-5.*(range+masl)+1.71*10^-9.*(range+masl).^2; 
%% Eta
% signal=NaN(nrange,nbins,dwells);
% ze=NaN(nrange,nbins,dwells);
eta = NaN(nrange,nbins,dwells);
%filename = '20151229_0000.mmclxa'; npw1nc=ncread([dpath,filename],'npw1');
vlx = NaN(nrange,dwells); vvx = NaN(nrange,dwells);
for nd = 1:dwells
    noise = 20.78*npw1(nd)+1831;%HSDV_co(rg,nd)*nbins
    for rg = 1:nrange
        [snr, v, vv]=calc_snr(f_co(rg,:,nd),noise); 
        vlx(rg,nd)=interp_jairo(v,1:nbins,vel);
        vvx(rg,nd)=interp_jairo(vv,1:nbins,vel);
        eta(rg,:,nd) = snr*RC(nd)*COFA(rg,nd)*((range(rg))/5000)^2*10^-18*pi^5*kwq/lambda^4;         
    end
end
% %% for Ze
% 
Ze=NaN(nrange,dwells);
for nd = 1:dwells    
    for rg = 1:nrange
        Ze(rg,nd)=sum(eta(rg,:,nd),'omitnan')/(10^-18*pi^5*kwq/lambda^4);
    end
end
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
t_min=15-6/1000*range(end);
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

% 'dD',dD,
dropsize = struct('drops',D,'dsd',dsd,'vel',vel,'Ze',Ze,...
    'vlx',vlx,'vvx',vvx,'fix',fix);
return

