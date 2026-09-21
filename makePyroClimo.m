function [zLon,zLat,zRad,zVec,zConf,frpAnn,frpDay,frpFrac,hLeg] = makePyroClimo(longBnds,latiBnds,ftVec,fRad,fConf,fLon,fLat,years2show,yearColor,tstr,legEnt,legPos,useConf,wantFRP,endYr,wantNorm)
% function [zLon,zLat,zRad,zVec,zConf,frpAnn,frpDay,frpFrac,hLeg] = makePyroClimo(longBnds,latiBnds,ftVec,fRad,fConf,fLon,fLat,years2show,yearColor,tstr,legEnt,legPos,useConf,wantFRP,endYr,wantNorm)
% Pyroclimographs (Hatchett et al. 2026) produced using satellite fire detections.
% Can show all detections using nominal confidence or user-selected as well as specific years with desired colors. 
%
% Input:
%   domain of interest (polygon defined by list of longBnds, latiBnds)
%   fLat,fLon,fRad,ftVec,fConf are latitude, longitude, fire radiative
%   power, and confidence of detection from satellite platform.
%   years2show: special years to highlight
%   yearColor: colors of years (in order)
%   tstr: title string
%   legEnt: legend entries (can state year, or year with specific fire
%   name)
%   legPos: legend position on the plot.
%   useConf: subset by n\% confidence (for full inclusion, NASA guidance
%       says use all detections)
%   wantFRP: 0 ignore, 1 if want include cumulative total FRP on right hand y-axis 
%   endYear: 0 ignore, 1 if year is complete but no detections, fill out with zeros to
%       end (last detect 1 Dec).
%   wantNorm: normalize counts by number of unique days
%
% Output:
%   zLon,zLat,zRad,zVec,zConf
% figure includes a two-panel plot with cumulative FRP by year (top) and
% total counts (bottom). Recommend a 37 unit wide by 22 unit high figure).
%
% Reference:
%   Hatchett, B.J., Lindley, T.T., Abatzoglou, J., Baring, A., Campbell, I., Gershunov, A., Guirguis, K., 
%       McCleod, G. Munroe, R., Nauslar, N.J., North, J.S., Orland, E., Rhoades, A.M. Short, K.C., Wells, E.M., 
%       Worsnop, R. P., 2026: Visualizing Pyroclimotology. Preprint submitted to EarthArXiv. https://doi.org/10.31223/X5JJ4D
%
%% first subset by region provided
[in,on]=inpolygon(fLon(:,1),fLat(:,1),longBnds,latiBnds);
in=in+on;
zLon = fLon(in==1);
zLat = fLat(in==1);
zVec = ftVec(in==1,:);
zConf = fConf(in==1);
zRad = fRad(in==1)./1000; % convert MW to GW;

if useConf

    foo = length(zRad);
    zRad = zRad(zConf>=30);
    zLon = zLon(zConf>=30);
    zLat = zLat(zConf>=30);
    zVec = zVec(zConf>=30,:);
    zConf = zConf(zConf>=30);
    dre = length(zRad);
    disp([num2str((dre/foo).*100) '% of values retained above nominal confidence detection']);
end

[a,b] = sort(zVec(:,8),'ascend');
frp = zRad(b,1);
frp(:,2) = zVec(b,8);
cumulFrp = cumsum(frp(:,1));
ind = unique(frp(:,2));
idx=find(~isnan(ind),1,'last');

idy = [];
for ii=1:idx
    idy(ii,1) = sum(frp(find(frp(:,2)==ind(ii)),1));
    idy(ii,2) = ind(ii);
end
[~,b] = sort(idy(:,1)','descend');
foo = [];
foo(:,1)=idy(b,1);
foo(:,2)=idy(b,2);
foo(:,3)=foo(:,1)./sum(foo(:,1));
foo(:,4)=cumsum(foo(:,1)./sum(foo(:,1)));

frpFrac = foo;
yrs = length(unique(zVec(:,1)));

clf
subplot(2,1,1)
ct = 1;
yrStart = min(zVec(:,1));
if endYr
    yrEnd = endYr;
else
    yrEnd = max(zVec(:,1));
end

frpAnn = zeros(26,1);
for iYear = yrStart:yrEnd
    if ismember(iYear,[2000 2004 2008 2012 2016 2020 2024])
        fulldVecs = datenum(2000,1,1):datenum(2000,12,31);
        fullVals = zeros(366,1);
        ind = find(zVec(:,7)>datenum(iYear,1,1,0,0,0) & zVec(:,7)<=datenum(iYear,12,31,23,59,0));
        dVec = [datenum(2000,zVec(ind,2),zVec(ind,3))];
    else
        fulldVecs = datenum(2001,1,1):datenum(2001,12,31);
        fullVals = zeros(365,1);
        ind = find(zVec(:,7)>datenum(iYear,1,1,0,0,0) & zVec(:,7)<=datenum(iYear,12,31,23,59,0));
        dVec = [datenum(2001,zVec(ind,2),zVec(ind,3))];
    end

    for ii=1:length(fullVals)
        if ismember(fulldVecs(ii),dVec)
            idx = find(dVec==fulldVecs(ii));
            fullVals(ii)=sum(zRad(ind(idx)));

        end
    end
    
    plotVals = zeros(365,1);

    % to plot leap years, distribute half of Feb 29 across day before/after
    % of 'normal' years. will need to update.
    if ismember(iYear,[2000 2004 2008 2012 2016 2020 2024 2028])
        plotVals(1:59) = fullVals(1:59);
        plotVals(60:365) = fullVals(61:366);
        plotVals(59) = fullVals(59)+fullVals(60).*0.5;
        plotVals(60) = fullVals(61)+fullVals(60).*0.5;
    else
        plotVals = fullVals;
    end

    fullVals = cumsum(plotVals);

    hold on;

    val = cumsum([zRad(ind)]);
    val(isnan(val))=0;
    if val
        frpAnn(iYear-(yrStart-1),1)=val(end);
    else % no val values!
        val = 0; % pad it
        frpAnn(iYear-(yrStart-1),1)=0;
    end
    
    fulldVecs = datenum(2001,1,1):datenum(2001,12,31);
    if iYear==2026
        try
            idy = find(dVec(end)<fulldVecs,1,'first')-1;
            plot(fulldVecs(1:idy),(fullVals(1:idy)),'linewi',2.5,'color',rgb('crimson'));
        catch
            idy = julday(2026,7,4); % fix at Jul 4 last upload
            plot(fulldVecs(1:idy),(fullVals(1:idy)),'linewi',2.5,'color',rgb('crimson'));
        end
        
        text(datenum(2002,1,21),(val(end)),num2str(iYear),'fontsize',13,'fontweight','bold','color',rgb('crimson'));

    else
        plot(fulldVecs,(fullVals),'linewi',1,'color',rgb('black'));
        text(datenum(2002,1,2),(val(end)),num2str(iYear),'fontsize',14)
    end
end

datetick('x','mmm')
xlim([datenum(2001,1,1) datenum(2001,12,31,12,0,0)])
ylabel({'Cumulative Fire Radiative'; 'Power [gigawatts]'},'FontWeight','bold');
xlabel('Month','FontWeight','bold')

grid on
dialFigR(17);
title(tstr)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Bar chart style timeseries by calendar day
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
hold on

% add julian days to zVec
zVec(:,9) = datenum(zVec(:,1),zVec(:,2),zVec(:,3))-datenum(zVec(:,1),1,1)+1;

for ii = 1:366
    ind = find(zVec(:,9)==ii);
    jday(ii,1) = sum(zRad(ind,1));
    jday(ii,2) = length(ind);
    jday(ii,3) = length(unique(zVec(ind,8)));
end

baz = 2; % baz = 1 for FRP; baz = 2 for counts
if wantNorm
    h=bar(1:366,jday(:,baz)./jday(ii,3),'facecolor',rgb('black'),'facealpha',0.2,'edgecolor','none');
else
    h=bar(1:366,jday(:,baz),'facecolor',rgb('black'),'facealpha',0.2,'edgecolor','none');
end
set(h,'barwidth',1.0);
frpDay = jday(:,2);
hold on

fday = NaN(366,24,2);
ct = 1; % count through days
ctt = 2; % count through colors, was set to one before I reversed the order to plot all years first (to avoid hazy colors from background all years)
for iYear = yrStart:yrEnd

    for ii = 1:366
        ind = find(zVec(:,1)==iYear & zVec(:,9)==ii);
        fday(ii,ct,1) = sum(zRad(ind,1));
        fday(ii,ct,2) = length(ind);
        fday(ii,ct,3) = length(unique(zVec(ind,8)));
    end

    baz = 2; % baz = 1 for FRP; baz = 2 for counts

    if ismember(iYear,years2show)
        %h=bar(1:366,squeeze(fday(:,ct,baz)),'facecolor',yearColor(ctt,:),'facealpha',0.75,'edgecolor','none');%yearColor(ctt,:),'edgealpha',0.75);
        %if wantNorm
        %h=bar(1:366,squeeze(fday(:,ct,baz)),'facecolor',yearColor(ctt,:),'facealpha',0.75,'edgecolor','none');%yearColor(ctt,:),'edgealpha',0.75);
        %else
        h=bar(1:366,squeeze(fday(:,ct,baz)),'facecolor',yearColor(ctt,:),'facealpha',0.75,'edgecolor','none');%yearColor(ctt,:),'edgealpha',0.75);
        %end
        
        set(h,'barwidth',1.0)
        ctt = ctt + 1;
    else
        %bar(1:366,squeeze(fday(:,ct,baz)),'facecolor',rgb('black'),'edgecolor','none','facealpha',0.5)
    end
    ct = ct + 1;
end

ylabel({'Total Daily Fire Detections'},'FontWeight','bold');

if wantFRP
    yyaxis('right')
    plot(1:366,cumsum(jday(:,1)),'k','linewi',0.92);
    ylabel({'Cumulative Fire Radiative';' Power [gigawatts]'},'FontWeight','bold'); % two panel
end

xlim([1 366]);
monDays = [1 32 60 91   121   152   182   213   244   274   305   335];
%monTitles = {'January';'February';'March';'April';'May';'June';'July';'August';'September';'October';'November';'December'};
monShortTitles = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
set(gca,'xtick',monDays,'xticklabel',monShortTitles);
dialFigR(17);
xlabel('Month','FontWeight','bold');
grid on
hLeg = legend;
hLeg.String=legEnt;
set(hLeg,'fontsize',14);

if wantFRP==1 && size(years2show,2)==1
    hLeg.NumColumns = 3;
    legPos=[0.22 0.475 0.1 0.025]; % top
elseif wantFRP==1 && size(years2show,2)>1
    hLeg.NumColumns = size(years2show,2);
    legPos=[0.2 0.475 0.1 0.025]; % top
elseif wantFRP==0 && ~isempty(legPos) && size(years2show,2)>1
    hLeg.NumColumns = size(years2show,1);
else
    legPos=[0.9472 0.35 0.01 0.01]; % along far right hand side
end
hLeg.Position=legPos;
yy=ylim;

% Temporal concentration (see Hatchett et al. 2026):
text(1,-1*yy(2)*.18,['50% of FRP in ' sprintf('%3.2f',100*(find(foo(:,4)>0.5,1,'first')/(yrs*365))) '% of days (' num2str(find(foo(:,4)>0.5,1,'first')) ' days)'] )
text(1,-1*yy(2)*.24,['75% of FRP in ' sprintf('%3.2f',100*(find(foo(:,4)>0.75,1,'first')/(yrs*365))) '% of days (' num2str(find(foo(:,4)>0.75,1,'first')) ' days)'] )
text(1,-1*yy(2)*.30,['90% of FRP in ' sprintf('%3.2f',100.*(find(foo(:,4)>0.9,1,'first')./(yrs*365))) '% of days (' num2str(find(foo(:,4)>0.9,1,'first')) ' days)'] )