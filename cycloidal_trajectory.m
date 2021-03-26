%% cycloidal_trajectory.m
%  It calculates a cycloidal trajectory given the position in joint
%  space and the initial and final time
%  
%  Output:
%       It returns a vector of cycloidal trajectory components
function [q,d_q,dd_q,ddd_q] = cycloidal_trajectory(qi, qf, ti, tf, Ts) 
    delta_T = tf-ti;
    delta_q = qf-qi;
    t = linspace(ti, tf, (tf-ti)/Ts);

    q = delta_q * ((t-ti)/delta_T - 1/(2*pi)*sin(2*pi*(t-ti)/delta_T)) + qi;
    d_q = delta_q/delta_T * (1 - cos((2*pi*(t-ti))/delta_T));
    dd_q = (2*pi*delta_q)/delta_T^2 * sin((2*pi*(t-ti))/delta_T);
    ddd_q = (4*pi^2*delta_q)/delta_T^3 * cos((2*pi*(t-ti))/delta_T);
end
