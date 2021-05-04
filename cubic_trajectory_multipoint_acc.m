function [q,dq,ddq,ts] = cubic_trajectory_multipoint_acc(qk,tk,dqi,dqf)
    if size(qk)~=size(tk)
        error("points and time values should have the same length");
    end
    n = size(qk,2)-1;

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

    A = zeros(6);
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

    c = zeros(n,1);
    for i = 1:n+1
        if i == 1
            c(i) = 6*((qk(2)-qk(1))/Tk(i)-dqi);
        elseif i == n+1
            c(i) = 6*(dqf-(qk(i)-qk(i-1))/Tk(i-1));
        else
            c(i) = 6*(((qk(i+1)-qk(i))/Tk(i))-((qk(i)-qk(i-1))/Tk(i-1)));
        end
    end
            
    ddq = thomas_algorithm(A, c);
    
    a0 = zeros(n,1);
    a1 = zeros(n,1);
    a2 = zeros(n,1);
    a3 = zeros(n,1);
    
    for k = 1:n
        a0(k) = qk(k);
        a1(k) = (qk(k+1)-qk(k))/Tk(k)-((ddq(k+1)+2*ddq(k))/6)*Tk(k);
        a2(k) = ddq(k)/2;
        a3(k) = (ddq(k+1)-ddq(k))/(6*Tk(k));
    end

    sample = 0.01;
    ts = tk(1):sample:tk(n+1);

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


