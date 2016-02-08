function [I,h] = mdm_nii_read(nii_fn)
% function [I,h] = mdm_nii_read(nii_fn)

% Possibly unzip file
[nii_fn, tmp_path, tmp_fn] = mdm_nii_gunzip(nii_fn);
if (~exist(nii_fn, 'file')), error(['File does not exist: ' nii_fn]); end

% Load header
h = mdm_nii_read_header(nii_fn);

% Open file and move beyond header to read image data; handle extended header data
fid = fopen(nii_fn, 'r');
fseek(fid, h.sizeof_hdr, -1);
extension = fread(fid, 4, 'uint8');
if (extension(1) ~= 0)
    n_bytes = fread(fid, 1, 'uint32');
else
    n_bytes = 0;
end
eff_hdr_size = h.sizeof_hdr + 4 + n_bytes;

fseek(fid, eff_hdr_size, -1);

if (ftell(fid) ~= eff_hdr_size)
    fclose(fid);
    error('something went wrong'); 
end

% Get image volume dimensions
dim_x = h.dim(2);
dim_y = h.dim(3);
dim_z = h.dim(4);
n_dyn = max(1,h.dim(5));

% Image is RBG-uint8 if bitpix equals 24
if (h.bitpix == 24)
    if (h.dim(1) > 3), error('no more than 3 dimensions expected'); end
    I = fread(fid, 3*dim_x*dim_y*dim_z, '*uint8');
    I = reshape(I, 3, dim_x, dim_y, dim_z);
else
    I = fread(fid, dim_x*dim_y*dim_z*n_dyn, ...
        ['*' nifti_datatype(h.datatype)]);
    I = reshape(I, dim_x, dim_y, dim_z, n_dyn);
end
fclose(fid);

% Remove the temp file and its dir
if (~isempty(tmp_fn)), delete(tmp_fn); end
if (~isempty(tmp_path)), rmdir(tmp_path); end