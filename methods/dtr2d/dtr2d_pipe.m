function fn = dtr2d_pipe(s, paths, opt)
% function fn = dtr2d_pipe(s, paths, opt)
%
% s     - input structure
% paths - either a pathname or a path structure (see mdm_paths)
% opt   - (optional) options that drive the pipeline
%
% fn    - a cell array with filenames to generated nii files

fn = '';

if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);
opt = dtr2d_opt(opt);
paths = mdm_paths(paths);
msf_log(['Starting ' mfilename], opt);

% Prepare: mask etc
if (opt.do_mask)
    s = mdm_s_mask(s, @mio_mask_threshold, [], opt);
end

% Run the analysis
if (opt.do_data2fit)
    if opt.do_bootstrap
        msf_mkdir(fileparts(paths.ind_fn));
        ind = (opt.dtr2d.ind_start-0) + round(rand([s.xps.n-(opt.dtr2d.ind_start-1),1])*(s.xps.n-(opt.dtr2d.ind_start-0)));
        save(paths.ind_fn, 'ind');
%         ind_fn = mdm_ind_save(ind, paths.ind_fn);
%         load(paths.ind_fn)
        opt.bootstrap.ind = ind;
    end
    mdm_data2fit(@dtr2d_4d_data2fit, s, paths.mfs_fn, opt);
    if opt.do_shuffle
        nodes4d = dtr2d_4d_m2nodes(paths.mfs_fn_init, opt);
        mdm_data2fit_supp(@dtr2d_4d_initfit2regfit, s, paths.mfs_fn, opt, nodes4d);
    end
end

if (opt.do_datafit2chisq)   
    opt.dtd.ind_start = opt.dtr2d.ind_start;
    chisq_fn = mio_datafit2chisq(@dtr2d_1d_fit2data, s, paths.mfs_fn, paths.chisq_fn, opt);
end

if (opt.do_fit2param)
    mdm_fit2param(@dtr2d_4d_fit2param, paths.mfs_fn, paths.dps_fn, opt);
end

% Save nifti parameter maps    
if (opt.do_param2nii)
    fn = mdm_param2nii(paths.dps_fn, paths.nii_path, opt.dtr2d, opt);
end

% Save dtr2d pdf   
if (opt.do_m2pdf)
    dtr2d_m2pdf(paths, opt);
end



