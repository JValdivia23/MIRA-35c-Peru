function fig1 = show_mmclx(filename,jvarname,gpath)
% show_mmclx - Plot a variable from NetCDF file of MIRA 35C radar
%     
%     This function returns a graph with automatic fix on dimension and labels.
%     Introducing show_mmclx(file,jvarname,gpath):
%     The inputs file and jvarname, are the NetCDF and variable to plot,
%     gpath is 'graph path'. 
%     jvarname should be as is described in ncdisp command.
% 
% Created by: Jairo Valdivia  -IGP                               Feb - 2017


% sfig= 1; % To save figure set sfig = 1
% if sfig
%     if ~exist('gpath','var'), gpath = 'J:\Otros\Radar\Plots\'; end
% end
if ~exist('boonlib.m','file'), disp('boonlib library does not exist'); end

if strcmp(filename(end-2:end),'.gz')
    if ~exist(filename(1:end-3),'file')
        gunzip(filename)
    end
        filename = filename(1:end-3);
        variable = ncread(filename,jvarname);
        range = ncread(filename,'range');
        time = ncread(filename,'time');
        jdel=input(['Delete ',filename(end-19:end),' file?(0:no, 1:yes): ']);
        if jdel, delete(filename), end
else
    variable = ncread(filename,jvarname);
    range = ncread(filename,'range');
    time = ncread(filename,'time');
end
fig1 = plot_mira35c(time,range,variable,jvarname);

end