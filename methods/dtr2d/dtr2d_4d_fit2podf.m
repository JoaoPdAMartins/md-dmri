function odf = dtr2d_4d_fit2podf(mfs_fn, odf_fn, opt)

%% Prepare options
% opt = mdm_opt();
% opt = dtr2d_opt(opt);
% opt.dtr2d.bin_disomax = [threshold_d_iso_big*1e-9];
% opt.dtr2d.bin_disomin = [0];
% opt.dtr2d.bin_sddeltamax = [1];
% opt.dtr2d.bin_sddeltamin = [threshold_d_delta_thin^2];
Nbins = numel(opt.dtr2d.bin_disomax);
N_nodes = opt.dtr2d.odf_nnodes;
% opt.dtr2d.r2extrap = -55e-3;

%% Smooth ODF parameters
angles_path = fullfile(framework_directory,'tools','uvec','repulsion_angles_tri');
angles = load(fullfile(angles_path,num2str(N_nodes))); %250, 350, 500, 1000, 3994, or 15970
ODindex = .05; % Orientation dispersion index of the Watson distribution smoothing kernel
odf_smooth.n = N_nodes;
odf_smooth.x = sin(angles.theta).*cos(angles.phi);
odf_smooth.y = sin(angles.theta).*sin(angles.phi);
odf_smooth.z = cos(angles.theta);
odf_smooth.tri = angles.tri;
odf_smooth.kappa = 1/tan(ODindex*pi/2);
                    
%% Get basic info
mfs = mdm_mfs_load(fullfile(bootstrap_directory, '1', 'mfs.mat'));
m = mfs.m;
mask = mfs.mask;
[Nx, Ny, Nz, N_sol] = size(m);
ind = false(N_sol,1);
ind(2:nb_dimension+1:end) = 1;
nn = (N_sol-1)/(nb_dimension+1);

%% Pre-read all of bootstrap realizations
m_global = zeros([NBS, Nx, Ny, Nz, N_sol]);
parfor nBS = 1:NBS
    mfs = mdm_mfs_load(fullfile(fullfile(bootstrap_directory,num2str(nBS)), 'mfs.mat'));
    m_global(nBS,:,:,:,:) = mfs.m;
end

%% Prepare median ODFs
odf_bsmedian.n = N_nodes; %250, 350, 500, 1000, 3994, or 15970
odf_bsmedian.x = sin(angles.theta).*cos(angles.phi);
odf_bsmedian.y = sin(angles.theta).*sin(angles.phi);
odf_bsmedian.z = cos(angles.theta);
odf_bsmedian.tri = angles.tri;
odf_bsmedian.kappa = 1/tan(ODindex*pi/2);

w_bin = zeros(Nbins, Nx, Ny, Nz, N_nodes);
diso_bin = zeros(Nbins, Nx, Ny, Nz, N_nodes);
sqddelta_bin = zeros(Nbins, Nx, Ny, Nz, N_nodes);
r2_bin = zeros(Nbins, Nx, Ny, Nz, N_nodes);
t2_bin = zeros(Nbins, Nx, Ny, Nz, N_nodes);

%% Compute median ODFs
parfor vx = 1:Nx    
    partial_w_bin = zeros(Nbins, Ny, Nz, N_nodes);
    partial_diso_bin = zeros(Nbins, Ny, Nz, N_nodes);
    partial_sqddelta_bin = zeros(Nbins, Ny, Nz, N_nodes);
    partial_r2_bin = zeros(Nbins, Ny, Nz, N_nodes);
    partial_t2_bin = zeros(Nbins, Ny, Nz, N_nodes);
   
    for vy = 1:Ny
        for vz = 1:Nz
            if mask(vx,vy,vz)       
                odf_w_temp = zeros(N_nodes, Nbins, NBS);
                odf_diso_temp = zeros(N_nodes, Nbins, NBS);
                odf_sqddelta_temp = zeros(N_nodes, Nbins, NBS);
                odf_r2_temp = zeros(N_nodes, Nbins, NBS);
                odf_t2_temp = zeros(N_nodes, Nbins, NBS);
                
                for nBS = 1:NBS
                    m = squeeze(m_global(nBS,vx,vy,vz,:));
                    dpar = m(circshift(ind,0,1));
                    dperp = m(circshift(ind,1,1));
                    theta = m(circshift(ind,2,1));
                    phi = m(circshift(ind,3,1));
                    r2 = m(circshift(ind,4,1));
                    w = m(circshift(ind,nb_dimension,1));
                    if isfield(opt.dtr2d,'r2extrap') == 1
                        w = w.*exp(-r2*opt.dtr2d.r2extrap);
                    end
                    
                    diso = (dpar + 2*dperp)/3;
                    daniso = (dpar - dperp)/3;
                    ddelta = msf_notfinite2zero(daniso./diso);
                    sqdaniso = daniso.^2;
                    sqddelta = msf_notfinite2zero(sqdaniso./diso.^2);
                    dratio = msf_notfinite2zero(dpar./dperp);
                    [dxx,dyy,dzz,dxy,dxz,dyz] = dtr2d_pars2elements(dpar,dperp,theta,phi);
                    t2 = msf_notfinite2zero(1./r2);
                    
                    dtds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
                        'sqdaniso',sqdaniso,'sqddelta',sqddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r2',r2,'t2',t2);
                    
                    for nbin = 1:Nbins
                        ind_bin = false([nn 4]);
                        ind_bin(:,1) = diso >= opt.dtr2d.bin_disomin(nbin);
                        ind_bin(:,2) = diso <= opt.dtr2d.bin_disomax(nbin);
                        ind_bin(:,3) = sqddelta >= opt.dtr2d.bin_sddeltamin(nbin);
                        ind_bin(:,4) = sqddelta <= opt.dtr2d.bin_sddeltamax(nbin);
                        ind_bin = all(ind_bin,2);
                        
                        if nnz(ind_bin) == 0 % No bin solution found for that specific bootstrap realization
                            odf_w_temp(:,nbin,nBS) = NaN([1, N_nodes]);
                            odf_diso_temp(:,nbin,nBS) = NaN([1, N_nodes]);
                            odf_sqddelta_temp(:,nbin,nBS) = NaN([1, N_nodes]);
                            odf_r2_temp(:,nbin,nBS) = NaN([1, N_nodes]);
                            odf_t2_temp(:,nbin,nBS) = NaN([1, N_nodes]);
                        else
                            dtds_temp = dtds;
                            dtds_temp.w = dtds.w.*ind_bin;
                            
                            % Discrete ODF
                            odf_discrete = struct;
                            odf_discrete.x = sin(dtds_temp.theta).*cos(dtds_temp.phi);
                            odf_discrete.y = sin(dtds_temp.theta).*sin(dtds_temp.phi);
                            odf_discrete.z = cos(dtds_temp.theta);
                            odf_discrete.w = dtds_temp.w;
                            odf_discrete.diso = dtds_temp.diso;
                            odf_discrete.sqddelta = dtds_temp.sqddelta;
                            odf_discrete.r2 = dtds_temp.r2;
                            odf_discrete.t2 = dtds_temp.t2;
                            
                            % Projection of the discrete ODF onto the smooth grid
                            odf_vox = dist_podf_discrete2smooth(odf_discrete, odf_smooth);
                            
                            odf_w_temp(:,nbin,nBS) = odf_vox.w(:);
                            odf_diso_temp(:,nbin,nBS) = odf_vox.diso(:);
                            odf_sqddelta_temp(:,nbin,nBS) = odf_vox.sqddelta(:);
                            odf_r2_temp(:,nbin,nBS) = odf_vox.r2(:);
                            odf_t2_temp(:,nbin,nBS) = odf_vox.t2(:);
                        end
                    end  
                end
                
                for nbin = 1:Nbins
                    partial_w_bin(nbin,vy,vz,:) = nanmedian(squeeze(odf_w_temp(:,nbin,:)),2); % Discards NaNs before computation
                    partial_diso_bin(nbin,vy,vz,:) = nanmedian(squeeze(odf_diso_temp(:,nbin,:)),2);
                    partial_sqddelta_bin(nbin,vy,vz,:) = nanmedian(squeeze(odf_sqddelta_temp(:,nbin,:)),2);
                    partial_r2_bin(nbin,vy,vz,:) = nanmedian(squeeze(odf_r2_temp(:,nbin,:)),2);
                    partial_t2_bin(nbin,vy,vz,:) = nanmedian(squeeze(odf_t2_temp(:,nbin,:)),2);
                end
                
            end
        end
    end
    
    w_bin(:,vx,:,:,:) = partial_w_bin;
    diso_bin(:,vx,:,:,:) = partial_diso_bin;
    sqddelta_bin(:,vx,:,:,:) = partial_sqddelta_bin;
    r2_bin(:,vx,:,:,:) = partial_r2_bin;
    t2_bin(:,vx,:,:,:) = partial_t2_bin;
end

for nbin = 1:Nbins
    odf_bsmedian.w_bin{nbin} = squeeze(w_bin(nbin,:,:,:,:));
    odf_bsmedian.diso_bin{nbin} = squeeze(diso_bin(nbin,:,:,:,:));
    odf_bsmedian.sqddelta_bin{nbin} = squeeze(sqddelta_bin(nbin,:,:,:,:));
    odf_bsmedian.r2_bin{nbin} = squeeze(r2_bin(nbin,:,:,:,:));
    odf_bsmedian.t2_bin{nbin} = squeeze(t2_bin(nbin,:,:,:,:));
end

odf = odf_bsmedian;