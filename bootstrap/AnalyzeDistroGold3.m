function [fitresult, gof] = AnalyzeDistroGold3(data)

[fdata,xdata] = ecdf(data);
logxdata = log(xdata);
mu1 = log(13.4);
mu2 = log(30.5);
sigma1 = log(1.51);
sigma2 = log(1.26);



[xData, yData] = prepareCurveData( logxdata, fdata );

% Set up fittype and options.
ft = fittype( 'a*normcdf(x,mu1,sigma1)+(1-a)*normcdf(x,mu2,sigma2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 mu1 mu2 sigma1 sigma2];
opts.StartPoint = [0.6 mu1 mu2 sigma1 sigma2];
opts.Upper = [1 mu1 mu2 sigma1 sigma2];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'Distrofit' );
% h = plot( fitresult, xData, yData );
% legend( h, 'data', 'fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'logx', 'Interpreter', 'none' );
% ylabel( 'fdata', 'Interpreter', 'none' );
% grid on