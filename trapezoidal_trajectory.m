%% trapezoidal_trajectory.m
%  It calculates a trapezoidal trajectory given all the values that are
%  assumed to be already checked
function [q,d_q,dd_q] = trapezoidal_trajectory(ti, tf, qi, qf, tc, d_qc, Ts)
    % Acceleration phase
    ta = linspace(ti, tc, (tc-ti)/Ts);
    qa = qi+(d_qc/(2*tc))*ta.^2;
    d_qa = (d_qc/tc)*ta;
    dd_qa = (d_qc/tc)+zeros(1,length(ta));

    % Linear phase
    tl = linspace(tc, tf-tc, (tf-tc-tc)/Ts);
    ql = qi + d_qc*(tl-tc/2);
    d_ql= d_qc+zeros(1,length(tl));
    dd_ql = 0+zeros(1,length(tl));

    % Deceleration phase
    td = linspace(tf-tc, tf, tc/Ts);
    qd = qf-((d_qc/(2*tc))*(tf-td).^2);
    d_qd = -(d_qc*(2*td - 2*tf))/(2*tc);
    dd_qd = -d_qc/tc+zeros(1,length(td));

    q = [qa ql qd];
    d_q = [d_qa d_ql d_qd];
    dd_q = [dd_qa dd_ql dd_qd];
end