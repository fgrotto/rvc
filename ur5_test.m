qi = [0.2 -0.5 0.3 -2 2.5 3]'; 
qf = [-0.2 -1.9 0.3 1 1.5 -1.2]'; 
d_qi_d = 0;
d_qf_d = 0;
ti = 0;
Ts = 0.001;
delta_T = NaN; %5
alfa = 1/3;
beta = 1/5;

vd_max = 2;  vd_min = -vd_max;
ad_max = 40;  ad_min = -ad_max;
jd_max = 60;  jd_min = -jd_max;

TimeValues = ti:Ts:tf;
DimValues = 6;
DataPosition = [];
DataVelocity = [];

for i = 1:length(qi)
    [q,d_q,dd_q,ddd_q,tf] = double_s_trajectory(qi(i),qf(i),d_qi_d,d_qf_d,vd_max,vd_min,ad_max,ad_min,jd_max,jd_min,delta_T,alfa,beta,ti,Ts);
    DataPosition = [DataPosition; q];
    DataVelocity = [DataVelocity; d_q];
end

qd.time = TimeValues;
qd.signals.values = DataPosition';
qd.signals.dimensions = DimValues;

dotqd.time = TimeValues;
dotqd.signals.values = DataVelocity';
dotqd.signals.dimensions = DimValues;