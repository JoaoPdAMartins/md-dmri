function [axh,scatterh,contourh] = mplot_3d_scatterdens_plot(axh,points_d,contourpars,axpars)

lw_contour = 1;
lw_axes = 2;
lw_plot = .5;
fs_axes = 15;

x = points_d.x(:); 
y = points_d.y(:); 
z = points_d.z(:); 
a = points_d.a(:); 
% r = points_d.r(:);
% g = points_d.g(:);
% b = points_d.b(:);
% c = repmat(points_d.bright(:),[1 3]).*[r g b];

dist_s.x = linspace(axpars.xmin,axpars.xmax,contourpars.Nx)';
dist_s.y = linspace(axpars.ymin,axpars.ymax,contourpars.Ny)';
dist_s.z = linspace(axpars.ymin,axpars.ymax,contourpars.Nz)';
dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));
dist_s.zsigma = 2*abs(dist_s.z(2) - dist_s.z(1));

if ~isfield(axpars,'no_scatter')
%     scatterh = scatter3(x,y,z,a+eps,[0 0 0]);
    [~,scatterh] = mplot_densityscattercolor(x,y,z,{'markerSize',1;...
        'xsigma',dist_s.xsigma;'ysigma',dist_s.ysigma;'zsigma',dist_s.zsigma});
    set(scatterh,'LineWidth',lw_plot)
else
    scatterh = scatter3(0,0,0,eps,[1 1 1]);
end

hold on

points_d.n = numel(x);
points_d.x = x;
points_d.y = y;
points_d.w = a;

dist_s = dist_2d_discrete2smooth(points_d,dist_s);
clevels_norm = linspace(contourpars.minlevel,contourpars.maxlevel,contourpars.Nlevels);
clevels = max(dist_s.w(:))*clevels_norm;

C = contourc(dist_s.x,dist_s.y,dist_s.w',clevels);
clear dist_s
count = 1;
contourh = [];
while count < length(C)
    numxy = C(2,count);
    xtemp = C(1,count+(1:numxy));
    ytemp = C(2,count+(1:numxy));
    h = plot3(xtemp,ytemp,axpars.zmin*ones(size(xtemp)),'k-','LineWidth',lw_contour);
    contourh = [contourh; h];
    count = count + numxy + 1;
end

points_d.n = numel(x);
points_d.x = x;
points_d.y = z;
points_d.w = a;

dist_s.x = linspace(axpars.xmin,axpars.xmax,contourpars.Nx)';
dist_s.y = linspace(axpars.zmin,axpars.zmax,contourpars.Nz)';
dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));

dist_s = dist_2d_discrete2smooth(points_d,dist_s);
clevels_norm = linspace(contourpars.minlevel,contourpars.maxlevel,contourpars.Nlevels);
clevels = max(dist_s.w(:))*clevels_norm;

C = contourc(dist_s.x,dist_s.y,dist_s.w',clevels);
clear dist_s
count = 1;
while count < length(C)
    numxy = C(2,count);
    xtemp = C(1,count+(1:numxy));
    ytemp = C(2,count+(1:numxy));
    h = plot3(xtemp,axpars.ymax*ones(size(xtemp)),ytemp,'k-','LineWidth',lw_contour);
    contourh = [contourh; h];
    count = count + numxy + 1;
end

points_d.n = numel(x);
points_d.x = y;
points_d.y = z;
points_d.w = a;

dist_s.x = linspace(axpars.ymin,axpars.ymax,contourpars.Ny)';
dist_s.y = linspace(axpars.zmin,axpars.zmax,contourpars.Nz)';
dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));

dist_s = dist_2d_discrete2smooth(points_d,dist_s);
clevels_norm = linspace(contourpars.minlevel,contourpars.maxlevel,contourpars.Nlevels);
clevels = max(dist_s.w(:))*clevels_norm;

C = contourc(dist_s.x,dist_s.y,dist_s.w',clevels);
clear dist_s
count = 1;
while count < length(C)
    numxy = C(2,count);
    xtemp = C(1,count+(1:numxy));
    ytemp = C(2,count+(1:numxy));
    h = plot3(axpars.xmin*ones(size(xtemp)),xtemp,ytemp,'k-','LineWidth',lw_contour);
    contourh = [contourh; h];
    count = count + numxy + 1;
end

axis([axpars.xmin axpars.xmax axpars.ymin axpars.ymax axpars.zmin axpars.zmax])
view(30,30)
set(contourh,'Color',.5*[1 1 1])
set(axh,'LineWidth',lw_axes,'FontSize',fs_axes)
axis(axh,'square')