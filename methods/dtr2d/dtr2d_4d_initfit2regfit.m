function res = dtr2d_4d_initfit2regfit(s, mfs_fn, opt, nodes4d)
% function res = dtr2d_4d_initfit2regfit(s, mfs_fn, opt, nodes4d)

if (nargin < 3), opt = []; end

res = -1;

if opt.do_bootstrap
    ind = opt.bootstrap.ind;
else
    ind = opt.dtr2d.ind_start:s.xps.n;
end

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal, dtr2d_nodes) dtr2d_1d_initfit2regfit(signal, xps, opt, ind, dtr2d_nodes);
dummy = mio_fit_model_suppdata(f, s, mfs_fn, opt, nodes4d);

res = 1;