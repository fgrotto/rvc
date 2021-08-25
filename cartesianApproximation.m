function [Q, QD, QDD, T] = cartesianApproximation(Pg)
    tk = 1:5;
    W = [1 100 100 1 1];
    mu = 1;
    Ts = 0.001;

    Q = [];
    QD = [];
    QDD = [];
    T = [];

    for i=1:3
        [q,qd,qdd,t] = cubic_trajectory_multipoint_smooth(Pg(i,:),tk,mu,W, Ts);
        Q = [Q;q];
        QD = [QD;qd];
        QDD = [QDD;qdd];
        T = [T;t];
    end
end

