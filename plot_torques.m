function [] = plot_torques(tau,t_s,Ts,maxTau,name)
    ti = t_s(1);
    tf = size(tau,2)*Ts+ti;
    time = linspace(ti,tf,(tf-ti)/Ts);

    [r,~] = size(tau);

    figure;
    set(gcf,'Name',name,'NumberTitle', 'off');
    for i = 1:r
        subplot(r,1,i)
        plot(time,tau(i,:))
        hold on;
        xlabel('time (s)')
        title(['q', num2str(i)])
        yline(maxTau(i))
        yline(-maxTau(i))
    end
end