function fig1 = plot_mira35c(time,range,variable,jvarname)
% plot_mira35c - Plot a variable of Mira 35c radar
%     
%     This function returns a graph with automatic fix on dimension and labels.
%     Introducing plot_mira35c(time,range,variable,jvarname,gpath):
%     The inputs time, range, variable and jvarname, are 'x axis', 'y axis',
%     'variable to graph', 'name of variable' and 'graph path', respectively; 
%     jvarname should be as is described in ncdisp command.
% 
% Created by: Jairo Valdivia                                    Sep - 2016

sfig= 0; % To save figure set sfig = 1
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
switch jvarname
    case {'Z','Ze','Zg'}
%         strcmp(jvarname,'Z') || strcmp(jvarname,'Ze') || strcmp(jvarname,'Zg')
        variable = 10*log10(variable);
        paleta = jet;
        palimit = [-60 30];
        etiqy = 'dBZ [mm^6/m^3]';
        jtitle = [jvarname,' '];
    case {'VEL','VELg','VELcl','VELrain'}
%         strcmp(jvarname,'VEL') || strcmp(jvarname,'VELg') || strcmp(jvarname,'VELcl') || strcmp(jvarname,'VELrain') || strcmp(jvarname,'VELplank')
%         variable = variable;
        paleta =  boonlib('rgmap');
        palimit = [-11 11];
        etiqy = '[m/s]';
        jtitle = 'Radial Velocity ';
    case {'RMS','RMSg'}
%         strcmp(jvarname,'RMS') || strcmp(jvarname,'RMSg') || strcmp(jvarname,'RMSplank') || strcmp(jvarname,'RMSrain') || strcmp(jvarname,'RMScl')
%         variable = variable;
        paleta = boonlib('rbmap');
        palimit = [0 4];
        etiqy = 'RMS [m/s]';
        jtitle = 'Spectral Width ';
    case {'LDR','LDRg','LDRplank','LDRcl','LDRrain'}
%         strcmp(jvarname,'LDR') || strcmp(jvarname,'LDRg') || strcmp(jvarname,'LDRplank') || strcmp(jvarname,'LDRrain') || strcmp(jvarname,'LDRcl')
        variable = 10*log10(variable);
        paleta = jet;
        palimit = [-35 5];
        etiqy = 'LDR dB';
        jtitle = [jvarname,' '];
    case {'SNR','SNRg','SNRplank'}
%         strcmp(jvarname,'SNR') || strcmp(jvarname,'SNRg') || strcmp(jvarname,'SNRplank') || strcmp(jvarname,'SNRrain') || strcmp(jvarname,'SNRcl')
        variable = 10*log10(variable);
        paleta = jet;
        palimit = [-25 80];
        etiqy = 'SNR dB';
        jtitle = [jvarname,' '];
    case {'RR'}
        paleta = parula(10); pppal=NaN(341,3);
        rticks=[0.1,1,2,3,5,8,13,21];
        dy=[1,rticks*10+1,341,342];
        for i= 1:10
            pppal(dy(i):dy(i+1)-1,:)=repmat(paleta(i,:),[dy(i+1)-dy(i),1]);
        end
        paleta=[pppal]; %[0.15,0.15,0.15];
        palimit = [0 34];
        etiqy = 'RR mm/h';
        jtitle = [jvarname,' '];
        rrt=1;
    case {'LWC'}
        variable = log10(variable);
        paleta = jet;
        palimit = [0 3];
        etiqy = 'LWC log_{10}(g/m^3)';
        jtitle = [jvarname,' '];
end

fig1=figure;
pcolor(meshtime,meshrange,variable)  % Ordenar gráfica (x,y,z)
shading flat                            % Borde de grilla 'off'
colormap(paleta)                           % Tipo de paleta
cb=colorbar;                            % Mostrar paleta
if rrt, cb=colorbar('ticks',[rticks,34]); end
caxis(palimit)                         % Rango de paletas
ylim([249/1000 13])
ylabel('Height AGL [km]')                % Etiquetas: > y
ylabel(cb,etiqy)               % > Paleta
xlabel(['Universal Time (hours) - ',datestr(timeok(1),1)]);
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

% Guardar imagen
if sfig
print(fig1,'-r220',[gpath,[jvarname,'-'],num2str(datestr(timeok(1),'yymmdd_HH'))],'-dpng')
% (variable,'resolución','título','formato')
end
end
