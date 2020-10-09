function odf = dtr2d_1d_fit2podf(m, opt, odf_s)

if nargin < 3
    odf_s = dist_smooth_odf(opt);
end

Nbins = numel(opt.dtr2d.bin_disomax);

% Extract parameters from distribution
[dpar, dperp, theta, phi, r2, w] = dtr2d_1d_m2pars(m);
if isfield(opt.dtr2d,'r2extrap') == 1
    w = w .* exp(-r2 * opt.dtr2d.r2extrap);
end
%
t2      = msf_notfinite2zero(1 ./ r2);
diso    = (dpar + 2*dperp) / 3;
daniso  = (dpar - dperp) / 3;
ddelta  = msf_notfinite2zero(daniso ./ diso);
sdaniso = daniso.^2;
sddelta = msf_notfinite2zero(sdaniso ./ diso.^2);
dratio  = msf_notfinite2zero(dpar ./ dperp);
[dxx, dyy, dzz, dxy, dxz, dyz] = dtr2d_pars2elements(dpar, dperp, theta, phi);

dtds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
    'sdaniso',sdaniso,'sddelta',sddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r2',r2,'t2',t2);

% Compute ODF using a number of pre-defined bins
for nbin = 1:Nbins
    ind_bin = false([numel(w) 4]);
    ind_bin(:,1) = diso >= opt.dtr2d.bin_disomin(nbin);
    ind_bin(:,2) = diso <= opt.dtr2d.bin_disomax(nbin);
    ind_bin(:,3) = sddelta >= opt.dtr2d.bin_sddeltamin(nbin);
    ind_bin(:,4) = sddelta <= opt.dtr2d.bin_sddeltamax(nbin);
    ind_bin = all(ind_bin,2);
    
    if (nnz(ind_bin) == 0 || sum(dtds.w .* ind_bin) < opt.dtr2d.w_thresh * sum(dtds.w)) % No bin solution found for that specific bootstrap realization
        odf.bin{nbin}.w       = zeros(opt.dtr2d.odf_nnodes, 1);
        odf.bin{nbin}.diso    = NaN(opt.dtr2d.odf_nnodes, 1);
        odf.bin{nbin}.sddelta = NaN(opt.dtr2d.odf_nnodes, 1);
        odf.bin{nbin}.r2      = NaN(opt.dtr2d.odf_nnodes, 1);
        odf.bin{nbin}.t2      = NaN(opt.dtr2d.odf_nnodes, 1);
    else
%         dtds_temp   = dtds;
%         dtds_temp.w = dtds.w .* ind_bin;

        % Discrete ODF
        odf_discrete = dtds;
        %
        odf_discrete.x       = sin(dtds.theta) .* cos(dtds.phi);
        odf_discrete.y       = sin(dtds.theta) .* sin(dtds.phi);
        odf_discrete.z       = cos(dtds.theta);
        odf_discrete.w       = dtds.w .* ind_bin;
        odf_discrete.diso    = dtds.diso;
        odf_discrete.sddelta = dtds.sddelta;
        odf_discrete.r2      = dtds.r2;
        odf_discrete.t2      = dtds.t2;
        
        % Projection of the discrete ODF onto the smooth grid
        odf_vox = dist_podf_discrete2smooth(odf_discrete, odf_s);
        
        odf.bin{nbin}.w       = odf_vox.w(:);
        odf.bin{nbin}.pv      = sum(odf_discrete.w(:)) / sum(dtds.w(:));
        odf.bin{nbin}.diso    = odf_vox.diso(:);
        odf.bin{nbin}.sddelta = odf_vox.sddelta(:);
        odf.bin{nbin}.r2      = odf_vox.r2(:);
        odf.bin{nbin}.t2      = odf_vox.t2(:);
    end
end