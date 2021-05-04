close all; 
clear all;

qi_d = 0;
qf_d = 10;
d_qi_d = 0.1;
d_qf_d = 0.3;
ti = 1;
Ts = 0.01;
delta_T = NaN; %5
alfa = 1/3;
beta = 1/5;

vd_max = 3;  vd_min = -vd_max;
ad_max = 3;  ad_min = -ad_max;
jd_max = 7;  jd_min = -jd_max;

[q,d_q,dd_q,ddd_q,tf] = double_s_trajectory(qi_d,qf_d,d_qi_d,d_qf_d,vd_max,vd_min,ad_max,ad_min,jd_max,jd_min,delta_T,alfa,beta,ti,Ts);
plot_trajectory('Double S Trajectory', ti:Ts:tf, [q; d_q; dd_q; ddd_q]); % vd_max, ad_max, jd_max)