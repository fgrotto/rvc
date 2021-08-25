function [P, DP, DDP, DDDP, T, p, dp, ddp] = linearMotion(pi, pf, Ts, P, DP, DDP, DDDP, T)
    
    s = 0: Ts: norm(pf - pi);

    if size(T)== 0 
        s_tmp  = 0: Ts: norm(pf - pi);
    else
        s_tmp = T(length(T)) + Ts: Ts: T(length(T)) + norm(pf - pi);
    end
    
    p = pi + s.*((pf - pi)/(norm(pf - pi)));
    dp = ones(3, length(s)).*(pf - pi)/(norm(pf - pi));
    ddp = zeros(3, length(s));
    dddp = zeros(3, length(s));
    
    P = [P, p];
    DP = [DP, dp];
    DDP = [DDP, ddp];
    DDDP = [DDDP, dddp];
    T = [T, s_tmp];
    
end

