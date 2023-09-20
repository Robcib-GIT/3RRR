load('toPlot-GA.mat')
toPlotGA = toPlot;
load('toPlot-NN.mat')
toPlotNN = toPlot;
load('toPlot-MCD2.mat')
boxplot([toPlot(:,end) toPlotGA(:,end) toPlotNN(:,end)],...
    'Labels',{'Numerical Method', 'GA', 'NN'})
axis([0.5 3.5 -0.02 0.5]);
title('Average error for the three DK Models');
