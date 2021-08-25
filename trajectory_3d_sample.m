%% 3D Trajectory - linear and circular motion primitives

clear all;  clc;

p1 = [0 0 0]';
p2 = [1 0 0]';
p3 = [2 1 0]';
p4 = [2 1 2]';
p5 = [2 0 2]';

Pg = [p1, p2, p3, p4, p5];

Ts = 0.001;

% centre of the circular motion
c_1 = [1 1 0]';
c_2 = [2 1 1]';
C = [c_1, c_2];

% true and false parameters 
tag = [1,1];

% sequence of the primitives
motion = [1, 2, 2, 1];

[P, DP, DDP, DDDP, T] = pipeline(Pg, C, tag, motion, Ts);

% [P, DP, DDP, T] = linearMotion(p5, p4, Ts, P, DP, DDP, T);
% 
% [P, DP, DDP, T] = circularMotion(p4, p3, Ts, c_2, P, DP, DDP, T, true);
% 
% [P, DP, DDP, T] = circularMotion(p3, p2, Ts, c_1, P, DP, DDP, T, false);
% 
% [P, DP, DDP, T] = linearMotion(p2, p1, Ts, P, DP, DDP, T);
 
[Q, QD, QDD, Tc] = cartesianApproximation(Pg);
 
figure(1);


p = plot3(P(1,:), P(2,:), P(3,:));
p.LineWidth = 5;
title('Position');xlabel('x');ylabel('y');zlabel('z');
grid on
hold on
plot3(Pg(1,:), Pg(2,:), Pg(3,:),'o', 'MarkerSize',12);
plot3(C(1,:), C(2,:), C(3,:),'o', 'MarkerSize',12)
hold on;
p = plot3(Q(1,:), Q(2,:), Q(3,:));
p.LineWidth = 5;
legend("Trajectory", "Points", "Centers", " Approximated Trajectory")

figure(2);

p = plot3(DP(1,:), DP(2,:), DP(3,:));
p.LineWidth = 5;
title('Velocity');xlabel('x');ylabel('y');zlabel('z');
grid on
% hold on
%plot3(Pg(1,:), Pg(2,:), Pg(3,:),'o', 'MarkerSize',12);
hold on;
p = plot3(QD(1,:), QD(2,:), QD(3,:));
p.LineWidth = 5;
legend("Trajectory", " Approximated Trajectory")

figure(3);

p = plot3(DDP(1,:), DDP(2,:), DDP(3,:));
p.LineWidth = 5;
title('Acceleration');xlabel('x');ylabel('y');zlabel('z');
grid on
% hold on
%plot3(Pg(1,:), Pg(2,:), Pg(3,:),'o', 'MarkerSize',12);
hold on;
p = plot3(QDD(1,:), QDD(2,:), QDD(3,:));
p.LineWidth = 5;
legend("Trajectory", " Approximated Trajectory")

figure(4);

p = plot3(DDDP(1,:), DDDP(2,:), DDDP(3,:));
p.LineWidth = 5;
title('Jerk');xlabel('x');ylabel('y');zlabel('z');
grid on
hold on
%plot3(Pg(1,:), Pg(2,:), Pg(3,:),'o', 'MarkerSize',12);