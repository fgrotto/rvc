%% trapezoidal_with_vel_trajectory.m
%  It calculates a trapezoidal trajectory given all the values that are
%  assumed to be already checked
function [q,d_q,dd_q] = trapezoidal_with_vel_trajectory(ti, tf, qi, qf, d_qi, d_qf, tc, d_qc, Ts)
    % Acceleration phase
    ta = linspace(ti, tc, (tc-ti)/Ts);
    qa = qi+d_qi*(ta-ti)+(d_qc-d_qi)/(2*(tc-ti))*(ta-ti).^2;
    d_qa = d_qi+(d_qc-d_qi)/(tc-ti)*(ta-ti);
    dd_qa = (d_qc-d_qi)/(tc-ti)+zeros(1,length(ta));

    % Linear phase
    tl = linspace(tc, tf-tc, (tf-tc-tc)/Ts);
    ql = qi + d_qi*(tc-ti)/2+d_qc*(tl-ti-(tc-ti)/2).^2;
    d_ql= d_qc+zeros(1,length(tl));
    dd_ql = 0+zeros(1,length(tl));

    % Deceleration phase
    td = linspace(tf-tc, tf, tc/Ts);
    qd = qf-d_qf*(tf-td)-((d_qc-d_qf)/(2*(ti+tc)))*(tf-td).^2;
    d_qd = d_qf+((d_qc-d_qf)/(ti+tc))*(tf-td);
    dd_qd = -(d_qc-d_qf)/(ti+tc)+zeros(1,length(td));

    q = [qa ql qd];
    d_q = [d_qa d_ql d_qd];
    dd_q = [dd_qa dd_ql dd_qd];
end