function [fitresult, gof] = AnalyzeDistro(data)

[fdata,xdata] = ecdf(data);
logxdata = log(xdata);

[xData, yData] = prepareCurveData( logxdata, fdata );

% Set up fittype and options.
ft = fittype( 'a*normcdf(x,mu1,sigma1)+(1-a)*normcdf(x,mu2,sigma2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 2.7 0 0];
opts.StartPoint = [0.6 2 4 0.1 0.1];
opts.Upper = [1 3.2 5 10 10];

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