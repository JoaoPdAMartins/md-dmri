function median_odf = mdm_odf_median(bs_odf)
% function median_odf = mdm_odf_median(bs_odf)
%

ind_bs = find(~cellfun('isempty',bs_odf));
if numel(ind_bs) == 0
    warning('odf_bs is empty')
    return
end

sz = ones(1,4);
sz_temp = size(bs_odf{ind_bs(1)}.bin{1}.w);
sz(1:numel(sz_temp)) = sz_temp;

f = fieldnames(bs_odf{ind_bs(1)});

for c = 1:numel(f)
    for cbin = 1:numel(bs_odf{ind_bs(1)}.(f{c}))
        fbin = fieldnames(bs_odf{ind_bs(1)}.(f{c}){cbin});
        for cfbin = 1:numel(fbin)
            % Double-check that ODFs have less than 4 dims
            % as bootstraps are stored and averaged in the 5th dim
            if (size(bs_odf{ind_bs(1)}.(f{c}){cbin}.(fbin{cfbin}), 1) == sz(1) && ndims(bs_odf{ind_bs(1)}.(f{c}){cbin}.(fbin{cfbin}))<5)
                odftemp = zeros(sz(1),sz(2),sz(3),sz(4),numel(ind_bs));
                for nbs = 1:numel(ind_bs)
                    odftemp(:,:,:,:,nbs) = bs_odf{ind_bs(nbs)}.(f{c}){cbin}.(fbin{cfbin});
                end
                median_odf.(f{c}){cbin}.(fbin{cfbin}) = msf_notfinite2zero(nanmedian(odftemp,5));
            end
        end
    end
end