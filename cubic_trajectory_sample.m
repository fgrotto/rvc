clear all;
close all;

qi = 0;
qf = 3;
dot_qi = 0;
dot_qf = 0;
ti = 1;
tf = 2;
Ts = 0.01;
t = linspace(ti, tf, (tf-ti)/Ts);

[q,d_q,dd_q] = cubic_trajectory_siciliano(qi, dot_qi, qf, dot_qf, ti, tf, Ts);
plot_trajectory('Cubic polynomial trajectory (Siciliano)', t, [q; d_q; dd_q])

[q,d_q,dd_q] = cubic_trajectory(qi, dot_qi, qf, dot_qf, ti, tf, Ts);
plot_trajectory('Cubic polynomial trajectory', t, [q; d_q; dd_q])