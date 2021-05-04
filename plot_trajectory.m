%% plot_trajectory.m
%  It plots the given trajectories
%  Input:
%         - name name of the plot
%         - time vector of time istances
%         - trajectories matrix pxn where p is the number of trajectories (position, 
%             velocities, accelerations, jerk, snap) and n is the samples
%             in the corresponding trajectories (ex, [q; d_q; ddd_q])
function plot_trajectory(name, time, trajectories) % vd_max, ad_max, jd_max)
    [p,~] = size(trajectories);
    
    plot_number = 321;
    if p == 4
       plot_number = 411;
    end
    if p == 3
       plot_number = 311;
    end
    
    if p > 0
        figure('Name', name);
        subplot(plot_number);
        plot(time, trajectories(1,:));
        xlabel("time [s]");
        ylabel("position [rad]");
        title("Position Plot");
%         title(['Position Plot  vmax=',num2str(vd_max),'  amax=',num2str(ad_max),'  jmax=',num2str(jd_max)]);
        xlim([time(1) time(end)]);
    end
    
    if p > 1
        subplot(plot_number+1);
        plot(time, trajectories(2,:));
        xlabel("time [s]");
        ylabel("velocity [rad/s]");
        title("Veloctiy Plot");
        xlim([time(1) time(end)]);
    end
    
    if p > 2
        subplot(plot_number+2);
        plot(time, trajectories(3,:));
        xlabel("time [s]");
        ylabel("acceleration [rad/s^2]");
        title('Acceleration Plot');
        xlim([time(1) time(end)]);
    end
    
    if p > 3
        subplot(plot_number+3);
        plot(time, trajectories(4,:));
        xlabel("time [s]");
        ylabel("jerk [rad/s^3]");
        title('Jerk Plot');
        xlim([time(1) time(end)]);
    end
    
    if p > 4
        subplot(plot_number+4);
        plot(time, trajectories(5,:));
        xlabel("time [s]");
        ylabel("snap [rad/s^4]");
        title('Snap Plot');
        xlim([time(1) time(end)]);
    end
end