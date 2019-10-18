function [d, h] = mplot_densityscattercolor(x,y,z,param)

% default values
plotting = 1; colMap = jet(100); markerSize = 20;
xsigma = .01; ysigma = xsigma; zsigma = xsigma;
r = sqrt((range(x)/30)^2 + (range(y)/30)^2 + (range(z)/30)^2);

% rewrite default parameters if needed
if nargin == nargin(mfilename)
    for j = 1:size(param,1), eval([param{j,1},'= param{j,2};']); end
end

% initialize
[m,n] = size([x y z]); d = zeros(m,1);
waith = waitbar(0,'Coloring your points!');

for j = 1:m % loop over every point    
    d(j) = sum(1/(xsigma*sqrt(2*pi)).*exp(-(x-x(j)).^2/(2*xsigma^2)).*...
        1/(ysigma*sqrt(2*pi)).*exp(-(y-y(j)).^2/(2*ysigma^2)).*...
        1/(zsigma*sqrt(2*pi)).*exp(-(z-z(j)).^2/(2*zsigma^2)));
    
    waitbar(j/m);
end
area = pi*r^2;
d = d/area;
close(waith);

% normalize density measure
maxd = max(d); mind = min(d); d = (d-mind)/(maxd-mind);

if plotting % plotting option
    T = linspace(0,1,size(colMap,1));
    col = interp1(T, colMap, d);
    if     n == 3
        h = scatter3(x,y,z,markerSize,col,'filled');
    else
        fprintf('Cannot plot in %d dimensions. Input x, y, z coordinates of 3D space.\n',m)
    end
end

end