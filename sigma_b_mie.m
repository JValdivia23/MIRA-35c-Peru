% Calculo de sigma_b
D = 0:0.1:10; %mm
a = D./2;
lambda = 3e8/34.85e9;%m
k_0 = 2*pi/lambda;
x = k_0*a*1e-3;

n = 4.5251;
k = -2.5754;

q_b = NaN(size(D));
for id = 1:numel(D)
    result = mie(complex(n,k),x(id));
    q_b(id) = result(7);
end
sigma_b = q_b.*(pi*a.^2)./(4*pi);

m = complex(n,k);
km = (m^2-1)/(m^2+2);
sigma_b_R = pi*5/lambda^4*abs(km)^2*(D.*1e-3).^6;

figure(1)
loglog(D,sigma_b)
hold on
loglog(D,sigma_b_R)
hold off
xlabel('D [mm]')
grid on

