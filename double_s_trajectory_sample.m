close all; 
clear all;

qi_d = 0;
qf_d = -2;
d_qi_d = 0;
d_qf_d = 1;
ti = 5;
Ts = 0.001;

vd_max = 1;  vd_min = -vd_max;
ad_max = 40;  ad_min = -ad_max;
jd_max = 60;  jd_min = -jd_max;

[q,d_q,dd_q,ddd_q,tf] = double_s_trajectory(qi_d,qf_d,d_qi_d,d_qf_d,vd_max,vd_min,ad_max,ad_min,jd_max,jd_min,ti,Ts);
plot_trajectory('Double S Trajectory', ti:Ts:tf, [q; d_q; dd_q; ddd_q])