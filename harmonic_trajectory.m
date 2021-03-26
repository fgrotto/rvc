%% harmonic_trajectory.m
%  It calculates a harmonic trajectory given the position in joint
%  space and the initial and final time
%  
%  Output:
%       It returns a vector of harmonic trajectory components
function [q,d_q,dd_q,ddd_q] = harmonic_trajectory(qi, qf, ti, tf, Ts) 
    delta_T = tf-ti;
    delta_q = qf-qi;
    t = linspace(ti, tf, (tf-ti)/Ts);

    q = (delta_q/2)*(1 - cos(pi*(t - ti)/delta_T)) + qi;
    d_q = pi*delta_q/(2*delta_T) * sin(pi*(t - ti)/delta_T);
    dd_q = pi^2*delta_q/(2*delta_T^2) * cos(pi*(t - ti)/delta_T);
    ddd_q = -pi^3*delta_q/(2*delta_T^3) * sin(pi*(t - ti)/delta_T);
end
