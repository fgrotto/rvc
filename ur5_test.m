qi = [0.2 -0.5 0.3 -2 2.5 3]'; 
qf = [-0.2 -1.9 0.3 1 1.5 -1.2]'; 
ti = 0; 
tf = 10; 
Ts = 0.01; 

TimeValues = linspace(ti,tf,(tf - ti)/Ts);
DimValues = 6;
DataPosition = [];
DataVelocity = [];

for i = 1:length(qi)
    [p,v,~,~] = harmonic_trajectory(qi(i),qf(i),ti,tf,Ts);
    DataPosition = [DataPosition; p];
    DataVelocity = [DataVelocity; v];
end

qd.time = TimeValues;
qd.signals.values = DataPosition';
qd.signals.dimensions = DimValues;

dotqd.time = TimeValues;
dotqd.signals.values = DataVelocity';
dotqd.signals.dimensions = DimValues;