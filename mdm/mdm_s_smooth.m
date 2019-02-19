function s = mdm_s_smooth(s, filter_sigma, o_path, opt)
% function s = mdm_s_smooth(s, filter_sigma, o_path, opt)

if (all(ischar(s))), s = setfield([], 'nii_fn', s); end

% init
if (nargin < 2), filter_sigma = 0.4; end
if (nargin < 3), o_path = fileparts(s.nii_fn); end
if (nargin < 4), opt.present = 1; end
opt = mdm_opt(opt);
msf_log(['Starting ' mfilename], opt);    


% construct the filename
[~,name] = msf_fileparts(s.nii_fn);
nii_fn = fullfile(o_path, [name '_s' opt.nii_ext]);

if (exist(nii_fn, 'file') && (~opt.do_overwrite))
    disp(['Skipping, output file already exists: ' nii_fn]);
    s.nii_fn = nii_fn;
    return; 
end

% write the mask, don't care if we overwrite anything
[I,h]   = mdm_nii_read(s.nii_fn);
I       = mio_smooth_4d(single(I), filter_sigma, opt);

mdm_nii_write(single(I), nii_fn, h);

s.nii_fn = nii_fn;

% always overwrite the xps if we've written a new image
opt.do_overwrite = 1; 
mdm_xps_save(s.xps, mdm_xps_fn_from_nii_fn(s.nii_fn), opt);
