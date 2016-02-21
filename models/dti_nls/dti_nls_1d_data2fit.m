function m = dti_nls_1d_data2fit(signal, xps, opt)
% function m = dti_nls_1d_data2fit(signal, xps, opt)
%
% Yields a 1x7 vector 'm' with the fit parameters
% m(1)   - s0
% m(2:7) - cholesky decomposition of the diffusion tensor

signal = double(signal);

unit_to_SI = [max(signal) [1 1 1 1 1 1] * 1e-9];

    function m = t2m(t)
        
        m(1) = t(1); % s0
        
        C = [...
            t(2) t(3) t(4);
            0    t(5) t(6);
            0      0  t(7)];
        
        m(2:7) = dtd_3x3_to_1x6(C' * C);

        m = m .* unit_to_SI;
    end
        
    function s = my_1d_fit2data(t,varargin)
        m = t2m(t);
        s = dti_nls_1d_fit2data(m, xps);
    end

t_guess = [1 1 1 1  0  0  0];
t_lb    = [0 0 0 0 -9 -9 -9];
t_ub    = [2 9 9 9 +9 +9 +9];

t = lsqcurvefit(@my_1d_fit2data, t_guess, [], signal, ...
    t_lb, t_ub, opt.dti_nls.lsq_opts);

m = t2m(t);


if (opt.dti_nls.do_plot)
    signal_fit = dti_nls_1d_fit2data(m, xps);
    x = (1:numel(signal))';
    plot(x,signal,'.',x,signal_fit,'o');
    pause(0.05);
end

end