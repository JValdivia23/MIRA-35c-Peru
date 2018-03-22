% CALCULAR LOS INDICES COMPLEJOS DE REFRACCION
%
% De aquí solo hay que extraer la parte correspondiente a 34.85GHz.
% Y sólo usa la parte para agua.
% 
% Considera e_r = complex(e_1,e_2);
% m2 = e_r; % este es el indice the refracción al cuadrado
% m = sqrt(m2);
% n = real(m);
% k = imag(m);
% 
% Creo que con éste código te ahorras la interpolación.
% Saludos,
% 
% Danny
% 
% ===========================
% close all
% clear all
% clc
 
%freq = [3 10 35]*1e9;
freq = 34.85e9;
temp1 = -20:35;% -20 -10 0];
 
 
% Part 1 and 2
 
for ifreq = 1:numel(freq)
    l = 299792458/freq(ifreq)*1e2; % in cm.
    for it = 1:numel(temp1)
        for im = 1%:2
            t = temp1(im,it);
            switch(im)
                case (1), ... % Water
                    mat = 'Water';
                    e_s = 78.54*(1.0-4.579e-3*(t-25)+1.19e-5*(t-25)^2-2.8e-8*(t-25)^3);
                    e_inf = 5.27137 + 0.021647*t-0.00131198*t^2;
                    a = -16.8129/(t+273)+0.0609265;
                    l_s = 0.00033836*exp(2513.98/(t+273));
                    s = 12.5664e8;
                case (2), ... % Ice
                    mat = 'Ice';
                    e_s = 203.168+2.5*t+0.15*t^2;
                    e_inf = 3.168;
                    a = 0.288+0.0052*t+0.00023*t^2;
                    l_s = 0.0009990288*exp(13200/((t+273)*1.9869));
                    s = 1.26*exp(-12500/((t+273)*1.9869));
            end
            e_1= e_inf + ((e_s-e_inf)*(1+(l_s/l)^(1-a)*sin(a*pi/2)))/ ...
                (1+2*(l_s/l)^(1-a)*sin(a*pi/2)+(l_s/l)^(2*(1-a)));
    
            e_2=((e_s-e_inf)*(l_s/l)^(1-a)*cos(a*pi/2))/ ...
                (1+2*(l_s/l)^(1-a)*sin(a*pi/2)+(l_s/l)^(2*(1-a)))+s*l/18.8496e10;
            m2 = complex(e_1,e_2); m =sqrt(m2); n = real(m); k = imag(m);
            %disp (['e_r = (' num2str(e_1) ',' num2str(e_2) '), freq = ' ...
            %    num2str(freq(ifreq)/1e9) 'GHz, temp = ' num2str(t) 'C, ' mat]);
            %disp([m2, m, n, k]);
            val_n(it) = n;
            val_k(it) = k;
            val_kwq(it) = (abs((m2-1)/(m2+2)))^2;           
        end
    end
end

