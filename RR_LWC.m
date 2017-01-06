%   RR
vel=-dropsize.vel;vel(vel<0)=NaN;
dwells=length(chunk.UTC);
nbins = chunk.Process_Param.sft;
nrange = length(chunk.range);
dD = dropsize.dD;
Z = NaN(nrange,dwells);
RR = NaN(nrange,dwells);
LWC = NaN(nrange,dwells);
% dD = NaN(range,nbins);
% for rg = 1:nrange
%     for n = 1:nbins
%         dD(rg,n) = 0.16/(6.18*exp(-0.6.*D(rg,n))*nue(rg));
% %             if n== 89, dD(n) = D(rg,n)-0.1; else dD(n) = D(rg,n)-D(rg,n+1); end
%     end
% end
%%
for nd = 1:dwells
    for rg = 1:nrange
        Z(rg,nd)=sum(dropsize.drops(rg,:).^6.*dropsize.dsd(rg,:,nd).*dD(rg,:),'omitnan');
        RR(rg,nd)=sum(0.0006*pi*dropsize.drops(rg,:).^3.*-vel.*dropsize.dsd(rg,:,nd).*dD(rg,:),'omitnan');
        LWC(rg,nd)=sum(pi/6*dropsize.drops(rg,:).^3.*dropsize.dsd(rg,:,nd).*dD(rg,:),'omitnan');
    end
end


% %%
% % para probar la diferencia de los delta del diametro dD
% for n = 1:nbins-1
%     dD1(n) = 0.16/(6.18*exp(-0.6.*D(rg,n))*nue(rg));
%     if n== 89, dD(n) = D(rg,n)-0.1; else dD(n) = D(rg,n)-D(rg,n+1); end
% end
% figure
% plot(D(1,1:end-1),dD1)
% hold on, plot(D(1,:),dD), hold off
% figure; plot(D(2,:),dv_dd(2,:))
% hold on; plot(D(2,:),0.16./dD)