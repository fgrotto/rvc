%% cubic_trajectory.m
%  It calculates a cubic spline given the position and velocity in joint
%  space and the initial and final time
%  
%  Output:
%       It returns a vector with cubic trajectories
function [q,d_q,dd_q] = cubic_trajectory(qi, dot_qi, qf, dot_qf, ti, tf, Ts) 
   a0 = qi;
   a1 = dot_qi;
   a2 = (3*(qf-qi)-(2*dot_qi+dot_qf)*(tf-ti))/(tf-ti)^2;
   a3 = (-2*(qf-qi)+(dot_qi+dot_qf)*(tf-ti))/(tf-ti)^3;
   
   x = [a3; a2; a1; a0];
   
   t = linspace(ti, tf, (tf-ti)/Ts);
   
   t = t-ti;
   q = x(1)*t.^3 + x(2)*t.^2 + x(3)*t + x(4);
   d_q = 3*x(1)*t.^2 + 2*x(2)*t + x(3);
   dd_q = 6*x(1)*t + 2*x(2);
end
