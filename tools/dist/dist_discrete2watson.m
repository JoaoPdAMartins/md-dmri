function dist_s = dist_discrete2watson( dist_d, dist_s)
% function odf_s = dist_odf_discrete2smooth(odf_d,odf_s)
%
% Convert from discrete to smooth Watson distribution. Gaussian convolution.
%
% odf_d discrete odf with fields
% odf_d.n number of discrete components
% odf_d.x vector of x-values
% odf_d.y vector of y-values
% odf_d.w vector of weights
% odf_s smooth odf  with input fields
% odf_s.x vector of x-values
% odf_s.y vector of y-values
% odf_s.kappa std for Gaussian convolution
% odf_s output smooth odf with additional fields
% odf_s.w matrix of weights
% odf_s.verts vertices

[K.x_s,K.x_d] = ndgrid(dist_s.x,dist_d.x);
[K.y_s,K.y_d] = ndgrid(dist_s.y,dist_d.y);
[K.z_s,K.z_d] = ndgrid(dist_s.z,dist_d.z);

k = exp(dist_s.kappa*(K.x_s.*K.x_d + K.y_s.*K.y_d + K.z_s.*K.z_d).^2);
k = k / sum(k(:));

clear K

dist_s.w = k*dist_d.w;
dist_s.verts = repmat(dist_s.w,[1 3]).*[dist_s.x dist_s.y dist_s.z];
