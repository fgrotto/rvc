%% cubic_trajectory_siciliano.m
%  It calculates a cubic spline given the position and velocity in joint
%  space and the initial and final time
%  
%  Output:
%       It returns a vector with cubic trajectories
function [q,d_q,dd_q] = cubic_trajectory_siciliano(qi, dot_qi, qf, dot_qf, ti, tf, Ts) 
   A = [ ti^3   ti^2 ti 1;
        3*ti^2  2*ti 1  0;
        tf^3    tf^2 tf 1;
        3*tf^2  2*tf 1  0];
   B = [qi; dot_qi; qf; dot_qf];
   % Solve the linear system AX = B
   x = linsolve(A,B);
   
   t = linspace(ti, tf, (tf-ti)/Ts);
   
   q = x(1)*t.^3 + x(2)*t.^2 + x(3)*t + x(4);
   d_q = 3*x(1)*t.^2 + 2*x(2)*t + x(3);
   dd_q = 6*x(1)*t + 2*x(2);
end
