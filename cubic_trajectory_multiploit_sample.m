clear all;
close all;

qk = [5 20 -25 6 -5 10];
tk = [0 2 8 12 15 20];

dqi = 0;
dqf = 0;
Ts = 0.01;

mu = 0.1;
wk = [0 1 5 10 1 1];

[q,dq,ddq,ts] = cubic_trajectory_multipoint_smooth(qk,tk,mu,wk,Ts);
% [q,dq,ddq,ts] = cubic_trajectory_multipoint_acc(qk,tk,dqi,dqf);

figure;
subplot(311);
hold on;
plot(ts,q);
plot(tk,qk,"r.");
xlabel("time(s)");
ylabel("q(rad)");
grid on;

subplot(312);
plot(ts,dq);
xlabel("time(s)");
ylabel("dq(rad/s)");
grid on;

subplot(313);
plot(ts,ddq);
xlabel("time(s)");
ylabel("ddq(rad/s^2)");
grid on;