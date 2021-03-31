function [] = subplotFD(file_name)
% Plotting Force Displacement Curve
output = fopen(file_name,'r');
formatSpec = '%f %f';
sizeA = [2 Inf];
plotData = fscanf(output,formatSpec,sizeA);
fclose(output);
plot3 = subplot(2,2,3);
cla(plot3)
plot(plotData(1,:),abs(plotData(2,:)),'Marker','*','MarkerSize',3,'LineWidth',1)
set(get(gca,'ylabel'),'String','Force','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')
set(get(gca,'xlabel'),'String','Displacement','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')
hold on
end