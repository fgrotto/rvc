function [q,dq,ddq,ts] = cubic_trajectory_multipoint_smooth(qk,tk,mu,wk, Ts)
    if size(qk)~=size(tk)
        error("points and time values should have the same length");
    end
    n = size(qk,2)-1;

    W = diag(wk);
    lambda = (1-mu)/6*mu;

    Tk = zeros(n,1);

    for i = 1:n
        Tk(i) = tk(i+1)-tk(i);
    end

    Ad = zeros(n+1,1);
    Ad(1) = 2*Tk(1);
    Ad(n+1) = 2*Tk(n-1);

    for i =2:n
        Ad(i) = 2*(Tk(i-1)+Tk(i));
    end

    A2 = zeros(n-1,1);
    for i = 1:n-1
        A2(i) = Tk(i);
    end

    A = zeros(n+1);
    for i = 1:n
        A(i,i) = Ad(i);
        if i == 1
            A(i,i+1) = A2(i);
        elseif i == n
            A(i,i-1) = A2(i-1);
        else
            A(i,i-1) = A2(i-1);
            A(i,i+1) =  A2(i);
        end
    end
    A(n+1,n+1) = 2*Tk(n);
    A(n+1-1,n+1) = Tk(n);
    A(n+1,n+1-1) = Tk(n);

    Cd = zeros(n+1,1);
    Cs = zeros(n,1);

    for k = 1:n+1
        if k==1
            Cd(1) = -6/Tk(1);
        elseif k == n+1
            Cd(n+1) = -6/Tk(n);
        else
            Cd(k) = -(6/Tk(k-1) + 6/Tk(k));
        end
    end

    for k = 1:n
        Cs(k) = 6/Tk(k);
    end

    C = diag(Cd) + diag(Cs,-1) + diag(Cs,1);

    s = inv(W + lambda*C'*inv(A)*C)*W*qk';
    dds = thomas_algorithm(A, C*s)';

    for k = 1:n
        a0(k) = s(k);
        a1(k) = (s(k+1)-s(k))/Tk(k)-((dds(k+1)+2*dds(k))/6)*Tk(k);
        a2(k) = dds(k)/2;
        a3(k) = (dds(k+1)-dds(k))/(6*Tk(k));
    end

    ts = tk(1):Ts:tk(n+1);

    j = 1;
    tj = 2;
    q = zeros(size(ts,2),1);
    dq = zeros(size(ts,2),1);
    ddq = zeros(size(ts,2),1);
    
    for i=1:size(ts,2)
        if ts(i) >tk(tj)
            j = j+1;
            tj = tj+1;
        end
        t = ts(i)-tk(j);
        q(i) = a3(j)*(t)^3+a2(j)*(t)^2+a1(j)*(t)+a0(j);
        dq(i) = 3*a3(j)*(t)^2+2*a2(j)*(t)+a1(j);
        ddq(i) = 6*a3(j)*(t)+2*a2(j);
    end
end


