%% plot_trajectory.m
%  It plots the given trajectories
%  Input:
%         - name name of the plot
%         - time vector of time istances
%         - trajectories matrix pxn where p is the number of trajectories (position, 
%             velocities, accelerations, jerk, snap) and n is the samples
%             in the corresponding trajectories (ex, [q; d_q; ddd_q])
function plot_trajectory(name, time, trajectories)
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
    end
    
    if p > 1
        subplot(plot_number+1);
        plot(time, trajectories(2,:));
        xlabel("time [s]");
        ylabel("velocity [rad/q]");
        title("Veloctiy Plot");
    end
    
    if p > 2
        subplot(plot_number+2);
        plot(time, trajectories(3,:));
        xlabel("time [s]");
        ylabel("acceleration [rad/q^2]");
        title('Acceleration Plot');
    end
    
    if p > 3
        subplot(plot_number+3);
        plot(time, trajectories(4,:));
        xlabel("time [s]");
        ylabel("jerk [rad/q^3]");
        title('Jerk Plot');
    end
    
    if p > 4
        subplot(plot_number+4);
        plot(time, trajectories(5,:));
        xlabel("time [s]");
        ylabel("snap [rad/q^4]");
        title('Snap Plot');
    end
end