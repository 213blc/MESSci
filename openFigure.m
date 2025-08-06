function openFigure(width,height,units)
%%
% openFigure.m pops open a figure with the specified dimensions and units. Note these
% dimensions do not appear to be preserved during export in all cases.
% Example: openFigure(35,23,'centimeters');
%%

figure('Units',units,'position',[0 0 width height],'paperpositionmode','auto')