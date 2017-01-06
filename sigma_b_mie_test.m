% Calculo de sigma_b
D = 0:0.1:10; %mm
a = D./2;
%lambda = 3e8/34.85e9;%m
xmt = 34.85*10^9; %frecuency of MIRA 35c 
c_light = 2.99792458*10^8;
lambda = c_light/xmt;
% lambda = 10e-2; %m
k_0 = 2*pi/lambda;
x = k_0*a*1e-3;

% m2 = 80.255+24.313i;
% 
% n = real(sqrt(m2));
% k = imag(sqrt(m2));
n = 5.45;
k = 2.90;

q_b = NaN(size(D));
for id = 1:numel(D)
    %result = mie(complex(n,k),x(id));
    [~,~,~,~,q_b(id),~] = bhmie(x(id),complex(n,k),1);
    %q_b(id) = result(7);
end
sigma_b = q_b.*pi^2.*(a.*1e-3).^2.*4; % faltaba transformar el radio a metros!!!

m = complex(n,k);
km = (m^2-1)/(m^2+2);
km2= 0.9152;%(abs(km))^2
sigma_b_R = (pi^5/lambda^4)*km2*(D.*1e-3).^6;

figure(9)
semilogy(D,sigma_b)
hold on
semilogy(D,sigma_b_R)
hold off
xlabel('D [mm]')
ylabel('Backscatter Cross-Section [m^2]')
grid on
axis([0,10,1e-10 1e3])


