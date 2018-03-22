function [Nc, pia, k] = get_pia(Na,SIGe,dR,strRg)

% Na = MIRA(1).dsd;
% SIGe = Mie.Extinction;
% dR = 31.1784; %m
% strRg = 4;

[nrange,nbins,ntime] = size(Na);
pia = ones(nrange,ntime);
k = NaN(nrange,ntime);
Nc = NaN(size(Na));
for nr = strRg:nrange
    Np = Na(nr,:,:).*permute(repmat(pia(nr-1,:),[128 1]),[3 1 2]);
    kp = nansum(repmat(SIGe(nr,:,:),[1 1 ntime]).*Np,2);
    Nc(nr,:,:) = -Np.*repmat(log(1-2.*kp.*dR)./(2*kp.*dR),[1 nbins 1]);
    k(nr,:) = squeeze(nansum(repmat(SIGe(nr,:,:),[1 1 ntime]).*Nc(nr,:,:),2));
    pia(nr,:) = pia(nr-1,:).*exp(2*k(nr,:)*dR);
    np10 = find(abs(pia(nr,:)) > 10 & pia(nr,:) < 0);
    pia(nr,np10) = pia(nr-1,np10);
end
