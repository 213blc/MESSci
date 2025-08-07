function warnByCWAnonlin(fileIn,binSteps,cbartitle,wantNum,outputName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%        National Weather Service Watches/Warnings/Information By CWA via IEM
%        Standard plot, no highlighted regions/CWAs/states
%        Non-linear colormap
% Example:
%  warnByCWAnonlin('frw_issuance.txt',[1 10 100 1000],{'Fire Warnings';'13 Jan 2006-3 Aug 2025'},0,'fire_warnings_nonlin_3Aug2025')
%  
%  fileIn = 'frw_issuance.txt';                             % input file
%  binSteps = [1 10 100 1000];                              % colormap binSteps must be larger than the largest number of observations 
%  cbartitle = {'Fire Warnings';'13 Jan 2006-3 Aug 2025'};  % title of your colorbar
%  wantNum = 0;                                             % want numbers plotted (warning, won't look good!)
%  outputName = 'frw_warnings_nonlin_3Aug2025';             % output filename
%
%  Citations:       NOAA Fire Weather Testbed, Baring, A., Hatchett, B.J., Hoekstra, S., McMeeking, L., 
%                   Thiem, K., Tolby, Z., Vickery, J., Wells, E.M., 2025: Fire Weather Testbed 
%                   Evaluations #002–004: An End-to-End Evaluation of NOAA's Emerging Wildland Fire Detection 
%                   and Warning Capabilities. NOAA Technical Memorandum OAR GSL-71, 
%                   https://doi.org/10.25923/4pqf-7g49
%                   
%                   Hatchett, B.J., 2025: MESSci: MATLAB tools for Earth
%                   System Science, [git hub code, activate DOI on zenodo]
%
%                   Hatchett, B.J., 2025: MESSci: MATLAB tools for Earth
%                   System Science, EarthArXiv, [add doi]
%
%                   Hatchett, B.J., 2025: National Weather Service Fire Warnings by Weather Forecast Office 
%                   [Data set]. Zenodo. https://doi.org/10.5281/zenodo.16749557
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% requires:
% Mapping Toolbox
% openFigure.m
% Maps: usastatehi.shp (states), GUM_adm0.shp (GU), wm713fm9130.shp (PR),
%       w_05mr24.shp (CWAs), re04oc12.shp (NWS regions)
%
%%
%%%% To dos:
% add background ocean/land
% simplify CWAs to avoid line corner effects on export? seems ok
% add numbers to all inset maps (currently done in AdobeIll)
% better nonlinear colormap?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Suppress polyshape warning
warning('off','MATLAB:polyshape:repairedBySimplify')

% Prepare Figure
openFigure(35,23,'centimeters');
ax = usamap('conus');
cwa = shaperead('/Users/benjamin.hatchett/Documents/Maps/NWS_CWA/w_05mr24.shp','UseGeoCoords',true);
reg = shaperead('/Users/benjamin.hatchett/Documents/Maps/NWS_regions/re04oc12.shp','UseGeoCoords',true);
states = shaperead('usastatehi.shp','UseGeoCoords',true);

% Import IEM output file
% rfw = importdata('/Users/benjamin.hatchett/Documents/Projects/IEM/frw_issuance.txt');
warn = importdata(fileIn);
warnCWA = warn.textdata(2:end,1);

% set up colormap
bins = binSteps;
clim([1 max(bins)]);
cmap = cbrewer('seq','YlOrRd',length(bins));
colormap(cmap(2:end,:));

ct = [];
colct = [];

for ii = 1:length(cwa)

    wfo = cwa(ii).FULLSTAID;

    % Matt B with the save! https://www.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string
    cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
    cell_array=warnCWA;
    string=wfo;
    logical_cells = cellfun(cellfind(wfo),cell_array);

    % total up instances of counts due to doubles in IEM data
    data = sum(warn.data(find(logical_cells==1),1));

    colbin = find(data>=bins,1,'last')+1;

    if isempty(colbin)
        colbin = 1;
    end

    if data>0
        llat = cwa(ii).Lat;
        llon = cwa(ii).Lon;
        fillm(llat,llon,'FaceColor',cmap(colbin,:),'LineWidth',0.75,'EdgeColor',rgb('silver'),'FaceAlpha',1);

        if wantNum
            a=polyshape(llon,llat);
            [bx,by]=centroid(a);
            textm(by,bx-0.45,num2str(data),"Color","White","fontsize",6);%,'fontweight','bold');
        end
    else % no data
        llat = cwa(ii).Lat;
        llon = cwa(ii).Lon;
        fillm(llat,llon,'FaceColor','none','LineWidth',0.75,'EdgeColor',rgb('silver'),'FaceAlpha',1);
    end
end

% the states
for ii = 1:51
    llat = states(ii).Lat;
    llon = states(ii).Lon;
    fillm(llat,llon,'FaceColor','none','LineWidth',0.75,'EdgeColor',rgb('black'));
end

% NWS regions
for ii = 1:6
    llat = reg(ii).Lat;
    llon = reg(ii).Lon;
    plotm(llat,llon,'LineWidth',2,'Color',rgb('black'))
end

gridm('off')
plabel off
mlabel off
gridm off

% Colorbar
cbar = colorbar('horizontal');
set(cbar,'position',[0.36 0.135 0.55 0.065],'xtick',0:(binSteps(end)/(length(binSteps)-1)):binSteps(end),'xticklabel',binSteps,'fontsize',16,'fontweight','bold','ticklength',0)
clim([0 binSteps(end)])

xlabel(cbar,cbartitle,'FontSize',18,'fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alaska
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2 = axes('Position',[0.08 0.15 0.28 0.28]);
usamap({'AK'})
setm(h2,'FFaceColor','w')
plabel off
mlabel off
gridm off
for ii = 123:125
    wfo = cwa(ii).FULLSTAID;

    % Matt B with the save! https://www.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string
    cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
    cell_array=warnCWA;
    string=wfo;
    logical_cells = cellfun(cellfind(wfo),cell_array);

    % total up instances of counts due to doubles in IEM data
    data = sum(warn.data(find(logical_cells==1),1));

    colbin = find(data>=bins,1,'last')+1;

    if isempty(colbin)
        colbin = 1;
    end

    if data>0
        llat = cwa(ii).Lat;
        llon = cwa(ii).Lon;
        fillm(llat,llon,'FaceColor',cmap(colbin,:),'LineWidth',0.75,'EdgeColor',rgb('silver'),'FaceAlpha',1);
        if wantNum
            a=polyshape(llon,llat);
            [bx,by]=centroid(a);
            textm(by,bx-0.85,num2str(data),"Color","White","fontsize",6);%,'fontweight','bold');
        end
    end
end

% regions
for ii = 1:6
    llat = reg(ii).Lat;
    llon = reg(ii).Lon;
    plotm(llat,llon,'LineWidth',2.2,'Color',rgb('black'))
end

% select CWAs
for ii = 1:length(ct)
    llat = cwa(ct(ii)).Lat;
    llon = cwa(ct(ii)).Lon;
    plotm(llat,llon,'LineWidth',3.1,'Color',rgb('black'))
end

% Hawaii
h2 = axes('Position',[0.75 0.46 0.12 0.12]);
axesm('MapProjection','mercator','MapLatLimit',[11 26],'MapLonLimit',[-160 -150])
setm(h2,'FFaceColor','w')
plabel off
mlabel off
gridm off
ii=46; % HFO
wfo = cwa(ii).FULLSTAID;

% Matt B with the save! https://www.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
cell_array=warnCWA;
string=wfo;
logical_cells = cellfun(cellfind(wfo),cell_array);

% total up instances of counts due to doubles in IEM data
data = sum(warn.data(find(logical_cells==1),1));

colbin = find(data>=bins,1,'last')+1;

if isempty(colbin)
    colbin = 1;
end

if data>0
    llat = cwa(ii).Lat;
    llon = cwa(ii).Lon;
    fillm(llat,llon,'FaceColor',cmap(colbin,:),'LineWidth',0.75,'EdgeColor',rgb('silver'),'FaceAlpha',1);
    if wantNum
        a=polyshape(llon,llat);
        [bx,by]=centroid(a);
        textm(by,bx-0.85,num2str(data),"Color","White","fontsize",6);%,'fontweight','bold');
    end
    
end

for ii = 1:6
    llat = reg(ii).Lat;
    llon = reg(ii).Lon;
    plotm(llat,llon,'LineWidth',2.2,'Color',rgb('black'))
end

for ii = 1:length(ct)
    llat = cwa(ct(ii)).Lat;
    llon = cwa(ct(ii)).Lon;
    fillm(llat,llon,'FaceColor',cmap(colct(ii),:),'LineWidth',0.5,'EdgeColor',rgb('black'),'FaceAlpha',1);
    plotm(llat,llon,'LineWidth',3.5,'Color',rgb('black'));
end

% Puerto Rico
h2 = axes('Position',[0.75 0.13 .125*0.8 0.25*0.8]);
axesm('MapProjection','mercator','MapLatLimit',[17.4 22],'MapLonLimit',[-68 -65])
setm(h2,'FFaceColor','w')
plabel off
mlabel off
gridm off
load coastlines

[m]=shaperead('/Users/benjamin.hatchett/Documents/Maps/puertorico/wm713fm9130.shp');
ii=105; % SJU
wfo = cwa(ii).FULLSTAID;

% Matt B with the save! https://www.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
cell_array=warnCWA;
string=wfo;
logical_cells = cellfun(cellfind(wfo),cell_array);

% total up instances of counts due to doubles in IEM data
data = sum(warn.data(find(logical_cells==1),1));

colbin = find(data>=bins,1,'last')+1;

if isempty(colbin)
    colbin = 1;
end

if data>0
    geoshow(cwa(ii),'FaceColor',cmap(colbin,:),'LineWidth',0.65,'EdgeColor',rgb('silver'),'FaceAlpha',0.9);
    if wantNum
        textm(19,-161,num2str(data),"Color","Black","fontsize",6);%,'fontweight','bold');
    end
end
geoshow(m.Y, m.X,'color','k','linewidth',2)

% Guam
h2 = axes('Position',[0.75 0.27 0.12*0.75 0.2*0.75]);
axesm('MapProjection','mercator','MapLatLimit',[12.9 13.8],'MapLonLimit',[144.4 145])
setm(h2,'FFaceColor','w')

[m]=shaperead('/Users/benjamin.hatchett/Documents/Maps/guam/GUM_adm0.shp');

ii=122; % Guam
wfo = cwa(ii).FULLSTAID;

% Matt B with the save! https://www.mathworks.com/matlabcentral/answers/2015-find-index-of-cells-containing-my-string
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
cell_array=warnCWA;
string=wfo;
logical_cells = cellfun(cellfind(wfo),cell_array);

% total up instances of counts due to doubles in IEM data
data = sum(warn.data(find(logical_cells==1),1));
colbin = find(data>=bins,1,'last')+1;

if isempty(colbin)
    colbin = 1;
end

if data>0
    geoshow(cwa(ii),'FaceColor',cmap(colbin,:),'LineWidth',0.65,'EdgeColor',rgb('silver'),'FaceAlpha',0.9);
    if wantNum
        textm(19,-161,num2str(data),"Color","Black","fontsize",6);%,'fontweight','bold');
    end
end

geoshow(m.Y, m.X,'color','k','linewidth',2)

print('-dpng','-vector','-r300',[outputName '.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Add point data, such as RAWS
% for ii = 1:length(raws)
%     if raws(ii,3)==2
%          plotm(raws(ii,1),raws(ii,2),'sw','MarkerFaceColor',rgb('dodgerblue'));
%     else
%         plotm(raws(ii,1),raws(ii,2),'pw','markerfacecolor','k','MarkerSize',10);
%     end
% end
%
% plotm(47.5,-78-4.5,'sw','markerfacecolor',rgb('dodgerblue'));
% textm(47.5,-76.5-4.5,'Shrub','fontsize',12)
%
% plotm(46.5,-78-4.7,'pw','markerfacecolor','k','MarkerSize',10);
% textm(46.5,-76.5-4.75,'Grassland','fontsize',12)
% %%