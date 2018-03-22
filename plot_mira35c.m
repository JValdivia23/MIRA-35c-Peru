function fig1 = plot_mira35c(time,range,variable,jvarname,gpath)
% plot_mira35c(time,range,variable,jvarname,gpath)
%
%     Plot a variable of Mira 35c radar
%     This function returns a graph with automatic fix on dimension and labels.
%     Introducing plot_mira35c(time,range,variable,jvarname,gpath):
%     The inputs time, range, variable and jvarname, are 'x axis', 'y axis',
%     'variable to graph', 'name of variable' and 'graph path', respectively; 
%     jvarname should be as is described in ncdisp command.
% 
% Created by: Jairo Valdivia                                    Sep - 2016

sfig= 0; % To autosave figure set sfig = 1
if exist('gpath','var'), sfig=1; end
if sfig
    if ~exist('gpath','var'), gpath = 'J:\Otros\Radar\Granizo\Graphs\'; end
end
if ~exist('boonlib.m','file'), disp('boonlib library does not exist'); end
%% 3.- Transformar el tiempo al lenguaje de fechas de matlab
% Matlab usa un vector cuya magnitud son días desde Enero, 00.00.0000 UTC
if time(1) > 1e6
    dayref=datenum([1970 01 01 00 00 00]);
    tdouble=double(time);
    timeok=tdouble/1440/60+dayref; % ! - UTC
else
    timeok=time;
end
%% 4.- Construir el arreglo de coordenadas
if range(end) > 1000
[meshtime, meshrange]=meshgrid(timeok,range/1000);
else
    [meshtime, meshrange]=meshgrid(timeok,range);
end
%% 5.- Graficar
rrt=0; % para ticks de RR
variable(1:2,:)=NaN; % Discarting noisy levels
switch jvarname
    case {'Z','Ze','Zg'}
        variable = 10*log10(variable);
        paleta = jet;
        palimit = [-40 60];
        etiqy = 'dBZ [mm^6/m^3]';
        jtitle = [jvarname,' '];
    case {'VEL','VELg','VELcl','VELrain'}
%         variable = variable;
        paleta =  boonlib('bjetmapx');%boonlib('rgmap');%
        palimit = [-10 10];
        etiqy = '[m/s]';
        jtitle = 'Radial Velocity ';
    case {'RMS','RMSg'}
%         variable = variable;
        paleta = boonlib('czmap',64);
        palimit = [-0.02 2.54];
        etiqy = 'RMS [m/s]';
        jtitle = 'Spectral Width ';
    case {'LDR','LDRg','LDRplank','LDRcl','LDRrain'}
        variable = 10*log10(variable);
        paleta = jet(16);
        palimit = [-35 5];
        etiqy = 'LDR dB';
        jtitle = [jvarname,' '];
    case {'SNR','SNRg','SNRplank'}
        variable = 10*log10(variable);
        paleta = jet;
        palimit = [-25 80];
        etiqy = 'SNR dB';
        jtitle = [jvarname,' '];
    case {'RR'}
        paleta = parula(10);
        rticks=[0,0.01,0.1,1,2,3,5,8,13,21,34];
        vari=NaN(size(variable));
        for x = 1:length(rticks)
            vari(variable>=rticks(x))=x;
        end
        variable=vari;
        paleta=[0.7,0.7,0.7; paleta];
        palimit = [1 12];
        etiqy = 'RR [mm/h]';
        jtitle = [jvarname,' '];
        rrt=1;
    case {'LWC'}
        variable = log10(variable);
        paleta = jet(64);
        palimit = [-1 3];
        etiqy = 'LWC log_{10}(mg/m^3)';
        jtitle = [jvarname,' '];
end

fig1=gcf;
pcolor(meshtime,meshrange,variable)  % Ordenar gráfica (x,y,z)
shading flat                            % Borde de grilla 'off'
colormap(paleta)                           % Tipo de paleta
cb=colorbar;                            % Mostrar paleta
if rrt, cb.TickLabels=[rticks(1:end-1),{'<34'},{''}]; end
caxis(palimit)                         % Rango de paletas
ylim([0 8])
ylabel('Height AGL [km]')                % Etiquetas: > y
ylabel(cb,etiqy)               % > Paleta
xlabel(['Time - ',datestr(timeok(1),1)]);
% xlabel('Universal Time')                    % > x
x1=(round(timeok(1)*24))/24;      % Indicar límite inferior
x2=(round(timeok(end)*24))/24;      % Indicar límite superior
xlim([x1 x2])                           % Cortar eje
% Rotar etiquetas
jblink=(x2-x1)*24; rotar90=0;
if jblink < 2
    jblink = 0.25;
elseif jblink >= 2 && jblink <4
    jblink = 0.5;
elseif jblink >= 4 && jblink <10
    jblink = 1;
elseif jblink >= 10
    jblink = 1; rotar90=1;
end
ax=gca; ax.XTick = x1:jblink/24:x2;          % Intervalos de una hora
if rotar90
    ax.XTickLabelRotation = 90;
end
dateaxis('x',15)                        % Formato de hora (ver 'help dateaxis' para detalles)
title([jtitle,datestr(timeok(1),1)]);                 % Título
fig1.CurrentAxes.FontSize=11;
fig1.CurrentAxes.LineWidth=1.5;
% Guardar imagen
if sfig
print(fig1,'-r300',[gpath,[jvarname,'-'],num2str(datestr(timeok(1),'yymmdd_HH'))],'-dpng')
disp(['Image saved in: ',gpath])
% (variable,'resolución','título','formato')
end
end
