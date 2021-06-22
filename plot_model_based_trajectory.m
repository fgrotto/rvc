function plot_model_based_trajectory(q,dq,ddq,Ts,q_s,t_s)
    [joints,~] = size(q);
    ti = t_s(1);
    tf = size(q,2)*Ts+ti;
    time = linspace(ti,tf,(tf-ti)/Ts);
    
    for i = 1:joints
       figure('NumberTitle', 'off', 'Name', ['Joint', num2str(i)]);
       subplot(311);
       plot(time,q(i,:));
       hold on;
       scatter(t_s, q_s, 'filled');
       xlabel("time [s]");
       ylabel("position [rad]");
       subplot(312);
       plot(time,dq(i,:));
       xlabel("time [s]");
       ylabel("velocity [rad/s]");
       subplot(313);
       plot(time,ddq(i,:));
       xlabel("time [s]");
       ylabel("acceleration [rad/s^2]");
       hold on;
    end
end
