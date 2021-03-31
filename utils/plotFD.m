function [] = plotFD(file_name)
% Plotting Force Displacement Curve

output = fopen(file_name,'r');
formatSpec = '%f %f';
sizeA = [2 Inf];
plotData = fscanf(output,formatSpec,sizeA);
fclose(output);
figure;
plot(plotData(1,:),abs(plotData(2,:)),'LineWidth',3)
set(get(gca,'ylabel'),'String','Force','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')
set(get(gca,'xlabel'),'String','Displacement','FontSize',16','FontWeight','bold','FontName','Times','Interpreter','tex')

end