function [t,ss] = msf_fit(f, s, l, u, n, opt)
% function [t,ss] = msf_fit(f, s, l, u, n, opt)
%
% Perform lsqcurvefit with function 'f' on data 's', repeating it 'n' times
%
% f - function handle
% s - signal
% l - lower bound
% u - upper bound
% n - number of repetitions
% opt - options for lsqcurvefit


ss = inf;
for c = 1:n
    
    g   = l + rand(size(l)) .* (u - l); % guess
    
    t_tmp = lsqcurvefit(f, g, [], s, l, u, opt); % tmp fit
    
    ss_tmp = sum( (s - f(t_tmp)).^2 );
    
    if (ss_tmp < ss)
        t = t_tmp;
        ss = ss_tmp;
    end
    
end
