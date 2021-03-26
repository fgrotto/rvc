%% fifth_order_trajectory_siciliano.m
%  It calculates a fifth order trajectory given the position,velocity and accelerations
%  in joint space and the initial and final time
%  
%  Output:
%       It returns a vector with fifth order trajectories
function [q,d_q,dd_q,ddd_q] = fifth_order_trajectory_siciliano(qi, dot_qi, dot_dot_qi, qf, dot_qf, dot_dot_qf, ti, tf, Ts) 
   A = [ ti^5   ti^4    ti^3    ti^2    ti  1;
        5*ti^4  4*ti^3  3*ti^2  2*ti    1   0;
        20*ti^3 12*ti^2 6*ti    2       0   0;
        tf^5    tf^4    tf^3    tf^2    tf  1;
        5*tf^4  4*tf^3  3*tf^2  2*tf    1   0;
        20*tf^3 12*tf^2 6*tf    2       0   0;];
   B = [qi; dot_qi; dot_dot_qi; qf; dot_qf; dot_dot_qf];
   % Solve the linear system AX = B
   x = linsolve(A,B);
   t = linspace(ti, tf, (tf-ti)/Ts);
   
   q = x(1)*t.^5 + x(2)*t.^4 + x(3)*t.^3 + x(4)*t.^2 + x(5)*t + x(6);
   d_q = 5*x(1)*t.^4 + 4*x(2)*t.^3 + 3*x(3)*t.^2 + 2*x(4)*t + x(5);
   dd_q = 20*x(1)*t.^3 + 12*x(2)*t.^2 + 6*x(3)*t + 2*x(4);
   ddd_q = 60*x(1)*t.^2+24*x(2)*t+6*x(3);
end

