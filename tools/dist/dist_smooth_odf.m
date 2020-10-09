function odf_s = dist_smooth_odf(opt)

run_path = mfilename('fullpath');
framework_path = fileparts(fileparts(fileparts(run_path)));
angles_path = fullfile(framework_path,'tools','uvec','repulsion_angles_tri');

% Possible choices: 250, 350, 500, 1000, 3994, or 15970
angles = load(fullfile(angles_path,num2str(opt.dtr2d.odf_nnodes)));

odf_s.n = opt.dtr2d.odf_nnodes;
odf_s.x = sin(angles.theta).*cos(angles.phi);
odf_s.y = sin(angles.theta).*sin(angles.phi);
odf_s.z = cos(angles.theta);
odf_s.tri = angles.tri;
odf_s.kappa = opt.dtr2d.odf_watsonkappa; % Orientation dispersion index of the Watson distribution smoothing kernel
% odf_s.kappa = 1/tan(ODindex*pi/2);