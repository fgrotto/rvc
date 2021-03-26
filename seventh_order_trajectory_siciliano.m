%% seventh_order_trajectory_siciliano.m
%  It calculates a fifth order trajectory given the position,velocity,accelerations
%  and jerks in joint space and the initial and final time
%  
%  Output:
%       It returns seventh order trajectories
function [q,d_q,dd_q,ddd_q,dddd_q] = seventh_order_trajectory_siciliano(qi, dqi, ddqi, dddqi, qf, dqf, ddqf, dddqf, ti, tf, Ts) 
   A = [ ti^7       ti^6        ti^5        ti^4    ti^3    ti^2    ti  1;
        7*ti^6      6^ti^5      5*ti^4      4*ti^3  3*ti^2  2*ti    1   0;
        42*ti^5     30*ti^4     20*ti^3     12*ti^2 6*ti    2       0   0;
        210*ti^4    120*ti^3    60*ti^2     24^ti   6       0       0   0;
        tf^7        tf^6         tf^5       tf^4    tf^3    tf^2    tf  1;
        7*tf^6      6^tf^5      5*tf^4      4*tf^3  3*tf^2  2*tf    1   0;
        42*tf^5     30*tf^4     20*tf^3     12*tf^2 6*tf    2       0   0;
        210*tf^4    120*tf^3    60*tf^2     24^tf   6       0       0   0;];
   B = [qi; dqi; ddqi;dddqi; qf; dqf; ddqf; dddqf];
   % Solve the linear system AX = B
   x = linsolve(A,B);
   t = linspace(ti, tf, (tf-ti)/Ts);
   
   q = x(1) *t.^7     +x(2) *  t.^6      +x(3) *  t.^5   +x(4) * t.^4  +x(5) *  t.^3   +x(6) * t.^2  +x(7) *  t +x(8) * 1;
   d_q = x(1)*7*t.^6   +x(2)*   6.^t.^5    +x(3) * 5*t.^4 +x(4) *4*t.^3 +x(5)* 3*t.^2 +x(6) *2*t  +x(7)*1;
   dd_q =  x(1)* 42*t.^5    +x(2)* 30*t.^4    +x(3)* 20*t.^3 +x(4)*12*t.^2 +x(5)*6*t   +x(6)* 2;
   ddd_q = x(1)*210*t.^4  +x(2)*120*t.^3 +  x(3)*60*t.^2 +x(4)*24.^t + x(5)* 6;
   dddd_q = 840*x(1)*t.^3 + 360*x(2)*t.^2 + 120*x(3)*t + 24*x(4);
end

