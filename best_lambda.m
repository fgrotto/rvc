function lambda = best_lambda(tau_m, tau_t)
    [r,c] = size(tau_m);
    if(r>c)
        tau_m = tau_m';
        r = c;
    end

    tmp = [];
    for i = 1:6
        tmp(i) = tau_m(i)/(max(abs(tau_t(i,:))));
    end

    lambda = min(tmp(:));
end