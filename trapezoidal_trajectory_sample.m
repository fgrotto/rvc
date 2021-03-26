%% Trapezoidal trajectory with null initial and final velocities
clear all;
close all;

Ts = 0.01;
ti = 0;
tf = 10;
qi = 0;
qf = 4;
tc = 1.2;
% d_qc = 0;
% dd_qc = 0;

% Case 1: Given tc(ti = 0)
dd_qc = (qf-qi)/((tc*tf)-tc^2);
d_qc = dd_qc*tc;

% Case 2: Given dd_qc(ti = 0)
% delta = (tf^2*dd_qc-4*(qf-qi))/dd_qc;
% if delta <= 0
%     error("delta must be > 0");
% end
% tc = tf/2-1/2*sqrt(delta);
% d_qc = dd_qc*tc;

% Case 3: Given d_qc(ti = 0)
% tc = (qi-qf+d_qc*tf)/d_qc;
% dd_qc = d_qc^2 / (qi-qf+d_qc*tf);
% if not(tc > 0)
%    error("tc must be greater than 0");
% end
% if not(tc<(tf-tc))
%    error("tc must be lower than tf-tc");
% end

% Case 4: Given dd_qc and d_qc (ti = 0)
% tc = d_qc/dd_qc;
% tf = (q_qc^2+dd_qc*(qf-qi))/(d_qc*dd_qc);
% if not((qf-qi) >= (d_qc^2/dd_qc))
%     error("a linear segment doesn't exist");
% end

t = linspace(ti, tf, (tf-ti)/Ts);
[q,d_q,dd_q] = trapezoidal_trajectory(ti,tf,qi,qf,tc,d_qc,Ts);

plot_trajectory('Trapezoidal trajectory', t, [q; d_q; dd_q]);

%% Trapezoidal trajectory with velocities
clear all; close all;

Ts = 0.01;
ti = 0;
tf = 10;
qi = 0;
qf = 2;
tc = 1;
d_qi = 0;
d_qf = 0;

dd_qc = (qf-qi)/((tc*tf)-tc^2);
d_qc = dd_qc*tc;

t = linspace(ti, tf, (tf-ti)/Ts);
[q,d_q,dd_q] = trapezoidal_with_vel_trajectory(ti, tf, qi, qf, d_qi, d_qf, tc, d_qc, Ts);

plot_trajectory('Trapezoidal trajectory', t, [q; d_q; dd_q]);