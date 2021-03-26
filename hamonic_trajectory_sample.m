clear all;
close all;

qi = 0;
qf = 10;
ti = 0;
tf = 8;
Ts = 0.01;

[q, d_q, dd_q, ddd_q] = harmonic_trajectory(qi,qf,ti,tf,Ts);

time = linspace(ti,tf,(tf-ti)/Ts);
plot_trajectory('Harmonic trajectory', time, [q; d_q; dd_q; ddd_q])


