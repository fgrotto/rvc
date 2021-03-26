clear all;
close all;

qi = 0;
qf = 2;
dot_qi = 0;
dot_qf = 0;
dot_dot_qi = 0;
dot_dot_qf = 0;
ddd_qi = 0;
ddd_qf = 0;
ti = 0;
tf = 8;
syms t;
Ts = 0.01;
t = linspace(ti, tf, (tf-ti)/Ts);

[q,d_q,dd_q,ddd_q,dddd_q] = seventh_order_trajectory_siciliano(qi, dot_qi, dot_dot_qi, ddd_qi, qf, dot_qf, dot_dot_qf, ddd_qf, ti, tf, Ts);

plot_trajectory('Seventh polynomial trajectory (Siciliano)', t, [q; d_q; dd_q; ddd_q; dddd_q])

[q,d_q,dd_q,ddd_q,dddd_q] = seventh_order_trajectory(qi, dot_qi, dot_dot_qi, ddd_qi, qf, dot_qf, dot_dot_qf, ddd_qf, ti, tf, Ts);

plot_trajectory('Seventh polynomial trajectory', t, [q; d_q; dd_q; ddd_q; dddd_q])