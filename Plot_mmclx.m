% =========================================================================
% La siguiente rutina esta orientada a una correcta gráfica de los archivos
% '.mmclx' del radar MIRA-35C
% CREADO POR: Jairo Valdivia - Setiembre 2016
% E-mail: valdiviaprado.ing@gmail.com
% =========================================================================

%% 1.- Rutas y nombre del archivo: (No olvidar el backslash '\' al final de la ruta)

dpath = 'J:\Otros\Radar\';  % RUTA DEL ARCHIVO
gpath = 'J:\Otros\Radar\Granizo\Graphs\'; % RUTA DE GRAFICOS
filename = '20160423_2200.mmclx'; % FILE
%%
Zg = [Zg ncread([dpath,filename],'Zg')];
VELg = [VELg ncread([dpath,filename],'VELg')];
time = [time; ncread([dpath,filename],'time')];
LDR = [LDR ncread([dpath,filename],'LDRg')];
RMS = [RMS ncread([dpath,filename],'RMSg')];
%% 2.- Llamar las variables de interés
% (Puedes ver las variables escribiendo "ncdisp([dpath,namefile])" en la
% ventana de comandos)

Zg = ncread([dpath,filename],'Zg');
%VELg = ncread([dpath,filename],'VELg');
range = ncread([dpath,filename],'range');
time = ncread([dpath,filename],'time');
%LDR = ncread([dpath,filename],'LDRg');
%RMS = ncread([dpath,filename],'RMSg');
% ! - el tiempo es un vector cuya magnitud son Segundos desde 01.01.1970 00:00 UTC

%% 3.- Transformar el tiempo al lenguaje de fechas de matlab
% Matlab usa un vector cuya magnitud son días desde Enero, 00.00.0000 UTC

dayref=datenum([1970 01 01 00 00 00]);
tdouble=double(time);
timeok=tdouble/1440/60+dayref; % ! - UTC

%% 4.- Construir el arreglo de coordenadas

[meshtime, meshrange]=meshgrid(timeok,range/1000);

%% 5.- Graficar

fig1=figure(1); 
pcolor(meshtime,meshrange,VELg)  % Ordenar gráfica (x,y,z)
shading flat                            % Borde de grilla 'off'
colormap(jet)                           % Tipo de paleta
cb=colorbar;                            % Mostrar paleta
caxis([-12 12])                         % Rango de paletas
ylim([range(3)/1000 13])
ylabel('Height AGL [km]')                % Etiquetas: > y
ylabel(cb,'dBZ [mm^6/m^3]')               % > Paleta
xlabel(['Universal Time (hours) - ',datestr(timeok(1),1)]);
% xlabel('Universal Time')                    % > x
x1=(floor(timeok(1)*24))/24;      % Indicar límite inferior
x2=(round(timeok(end)*24))/24;      % Indicar límite superior
xlim([x1 x2])                           % Cortar eje
% Rotar etiquetas
ax=gca; ax.XTick = x1:0.5/24:x2;          % Intervalos de una hora

dateaxis('x',15)                        % Formato de hora (ver 'help dateaxis' para detalles)
title(['Z ',datestr(timeok(1),1)]);                 % Título

% Guardar imagen
print(fig1,'-r0',[gpath,'Z_',num2str(datestr(timeok(1),'yymmdd_HH'))],'-dpng')
% (variable,'resolución','título','formato')

%% Arreglos de eje x 'tiempo'
% Solo para mejor visualización
% 
% % Limitar eje x
% x1=datenum([2016 01 06 2 00 00]);      % Indicar límite inferior
% x2=datenum([2016 01 06 7 00 00]);      % Indicar límite superior
% xlim([x1 x2])                           % Cortar eje
% % Rotar etiquetas
% ax=gca; ax.XTick = x1:1/24:x2;          % Intervalos de una hora
% ax.XTickLabelRotation = 90;             % Rotación de etiqueta
% % Formato
% dateaxis('x',15)                        % Formato de etiqueta
% 


