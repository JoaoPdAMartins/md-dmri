function odf_peaks = dtr2d_1d_peaks_from_podf(odf, opt, odf_s)

if nargin < 3
    odf_s = dist_smooth_odf(opt);
end

odf_s.w                = odf.w;
odf_s.w_normalized     = odf.w / sum(odf.w);
odf_s.diso             = odf.diso;
odf_s.sddelta          = odf.sddelta;
odf_s.r2               = odf.r2;
odf_s.t2               = odf.t2;
odf_s.verts            = repmat(odf_s.w,[1 3]).*[odf_s.x odf_s.y odf_s.z];
odf_s.verts_norm       = sqrt( sum( odf_s.verts.^2, 2));
odf_s.verts_normalized = bsxfun(@rdivide, odf_s.verts, odf_s.verts_norm);

% Define connectivity structure
TR = triangulation(odf_s.tri, odf_s.verts);
conn.basis = vertexAttachments(TR);
ind = false(size(conn.basis,1),1);
for i = 1:size(conn.basis,1)
    conn.tri = conn.basis{i};
    conn.verts = unique([odf_s.tri(conn.tri,1); odf_s.tri(conn.tri,2); odf_s.tri(conn.tri,3)]);
    if all(odf_s.verts_norm(i) >= odf_s.verts_norm(conn.verts))
        ind(i) = 1;
    end
end

% Filtering out low-probability peaks
indw = odf_s.w_normalized >= opt.dtr2d.peak_thresh * max(odf_s.w_normalized);
targ_f = {'verts', 'verts_normalized', 'diso', 'sddelta', 'r2', 't2', ...
    'w_normalized', 'w'};
for f = 1:numel(targ_f)
    odf_peaks.(targ_f{f}) = odf_s.(targ_f{f})(ind & indw, :);
end

% Filtering out redundant antipodal points (keep z > -.05)
ind_z_positive = odf_peaks.verts(:,3) > -.05;
for f = 1:numel(targ_f)
    odf_peaks.(targ_f{f}) = odf_peaks.(targ_f{f})(ind_z_positive, :);
end

% Filter out peaks with similar orientations
verts_inner = triu(odf_peaks.verts_normalized * odf_peaks.verts_normalized', 1);
ind_redund = abs(verts_inner) > abs( cos(opt.dtr2d.angle_thresh));
ind_redund = ~any(ind_redund, 1);
for f = 1:numel(targ_f)
    odf_peaks.(targ_f{f}) = odf_peaks.(targ_f{f})(ind_redund, :);
end

norm_odf_verts_peaks = sqrt(odf_peaks.verts(:,1).^2+odf_peaks.verts(:,2).^2 + odf_peaks.verts(:,3).^2);
[~, ind] = sort(norm_odf_verts_peaks,'descend');

% Select the n_peaks peaks with highest weights
odf_peaks.n = opt.dtr2d.n_peaks;
if numel(ind) < opt.dtr2d.n_peaks    
    for f = 1:numel(targ_f)
        odf_peaks.(targ_f{f}) = odf_peaks.(targ_f{f})(ind, :);
    end
    odf_peaks.norm(1:numel(ind),:) = norm_odf_verts_peaks(ind);
    odf_peaks.n = numel(ind);
else
    for f = 1:numel(targ_f)
        odf_peaks.(targ_f{f}) = odf_peaks.(targ_f{f})(ind(1:odf_peaks.n), :);
    end
    odf_peaks.norm = norm_odf_verts_peaks(ind(1:odf_peaks.n));    
end