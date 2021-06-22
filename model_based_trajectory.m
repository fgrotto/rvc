function [pos,vel,acc,tau,tau_raw, lambda,q_raw,dq_raw,ddq_raw] = model_based_trajectory(dh,q_s,t_s,maxTau,method,Ts)
    gravity = [0 0 9.81]';

    if (strcmp('euler', method))
        for i = 1:dh.dof
            [q,d_q,dd_q] = multipoint_trajectory(q_s, t_s, 0, 0, Ts);
            pos(i,:) = q;
            vel(i,:) = d_q;
            acc(i,:) = dd_q;
        end
    end

    if (strcmp('acc', method))
        for i = 1:dh.dof
            [q,d_q,dd_q] = multipoint_trajectory_acc_cont(q_s, t_s, 0, 0, Ts);

            pos(i,:) = q;
            vel(i,:) = d_q;
            acc(i,:) = dd_q;
        end
    end

    q_raw = pos;
    dq_raw = vel;
    ddq_raw = acc;

    tau_raw = [];

    for  i = 1:size(q,2)
        t = inv_dyn_recursive_NewtonEulero(dh, pos(:,i), vel(:,i), acc(:,i), gravity);
        tau_raw = [tau_raw,t];
    end

    lambda = sqrt(best_lambda(maxTau,tau_raw));
    t_s = t_s./lambda;
    
    pos = [];
    vel = [];
    acc = [];
    tau = [];
    
    if (strcmp('euler', method))
        for i = 1:dh.dof
            [q,d_q,dd_q] = multipoint_trajectory(q_s, t_s, 0, 0, Ts);
            pos(i,:) = q;
            vel(i,:) = d_q;
            acc(i,:) = dd_q;
        end
    end

    if (strcmp('acc', method))
        for i = 1:dh.dof
            [q,d_q,dd_q] = multipoint_trajectory_acc_cont(q_s, t_s, 0, 0, Ts);

            pos(i,:) = q;
            vel(i,:) = d_q;
            acc(i,:) = dd_q;
        end
    end
    
    for  i = 1:size(q,2)
        t = inv_dyn_recursive_NewtonEulero(dh, pos(:,i), vel(:,i), acc(:,i), gravity);
        tau = [tau,t];
    end
end
