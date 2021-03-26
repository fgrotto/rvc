clear all;
%close all;

qi=0;
qf=2;
d_qi=0;
d_qf=0;
d_qmax=1;
d_qmin=-d_qmax;
dd_qmax=2;
dd_qmin=-dd_qmax;
ddd_q_max=3;
ddd_q_min=ddd_q_max;
ti = 1;
Ts = 0.01;

[q,d_q,dd_q,ddd_q, tf]=double_s_trajectory(qi,qf,d_qi,d_qf,d_qmin,d_qmax,dd_qmin,dd_qmax,ddd_q_min,ddd_q_max,ti,Ts);
t = linspace(ti, tf, (tf-ti)/Ts);
plot_trajectory('Double S Trajectory', t, [q; d_q; dd_q; ddd_q])

function  [q,d_q,dd_q,ddd_q, tf]=double_s_trajectory(qi,qf,d_qi,d_qf,d_qmin,d_qmax,dd_qmin,dd_qmax,ddd_q_min,ddd_q_max, ti,Ts)
    sigma = sign(qf-qi);
    qi = sigma * qi;
    qf = sigma * qf;
    d_qi = sigma * d_qi;
    d_qf = sigma * d_qf;

    d_qmax = (sigma+1)/2*d_qmax + (sigma-1)/2*d_qmin;
%     d_qmin = (sigma+1)/2*d_qmin + (sigma-1)/2*d_qmax;
    dd_qmax = (sigma+1)/2*dd_qmax + (sigma-1)/2*dd_qmin;
%     dd_qmin = (sigma+1)/2*dd_qmin + (sigma-1)/2*dd_qmax;
    ddd_q_max = (sigma+1)/2*ddd_q_max + (sigma-1)/2*ddd_q_min;
%     ddd_q_min = (sigma+1)/2*ddd_q_min + (sigma-1)/2*ddd_q_max;

    % Case 1: d_qlim = d_qmax
    if (d_qmax-d_qi)*ddd_q_max<dd_qmax^2    
        % dd_qmax is not reached
        Tj1=sqrt((d_qmax-d_qi)/ddd_q_max);
        ta=2*Tj1;
    else                        
        % dd_qmax is actually reached
        Tj1=dd_qmax/ddd_q_max;
        ta=Tj1+(d_qmax-d_qi)/dd_qmax;
    end
    if (d_qmax-d_qf)*ddd_q_max<dd_qmax^2    
        % dd_qmin is not reached
        Tj2=sqrt((d_qmax-d_qf)/ddd_q_max);
        td=2*Tj2;
    else                       
        % dd_qmin is actually reached
        Tj2=dd_qmax/ddd_q_max;
        td=Tj2+(d_qmax-d_qf)/dd_qmax;
    end

    tv=(qf-qi)/d_qmax-ta/2*(1+d_qi/d_qmax)-td/2*(1+d_qf/d_qmax);

    if tv>0
        % The maximum velocity is actually reached
        d_qlim=d_qmax;
        dd_qlima=ddd_q_max*Tj1;
        dd_qlimd=-ddd_q_max*Tj2;
    end

    % Case 2 d_qlim < d_qmax
    if tv<=0 
        % This means that d_qlim < d_qmax
        tv=0;
        Tj1=dd_qmax/ddd_q_max;
        Tj2=Tj1;
        tj=Tj1;

        D=dd_qmax^4/(ddd_q_max^2)+2*(d_qi^2+d_qf^2)+dd_qmax*(4*(qf-qi)-2*(dd_qmax/ddd_q_max)*(d_qi+d_qf));
        ta=(dd_qmax^2/ddd_q_max-2*d_qi+sqrt(D))/(2*dd_qmax);
        td=(dd_qmax^2/ddd_q_max-2*d_qf+sqrt(D))/(2*dd_qmax);

        assert(ta<2*tj || td<2*tj, "the maximum (minimum) acceleration is not reached: decrease the value of dd_qmax");
        % assert(ta<0 || td<0, "only one of the acceleration or deceleration phase is necessary");
        if (ta<0 || td<0)
            display("only one of the acceleration or deceleration pahse is necessary")
        end

        dd_qlima=ddd_q_max*Tj1;
        dd_qlimd=-ddd_q_max*Tj2;
        d_qlim=d_qi+(ta-Tj1)*dd_qlima; % Same as d_qlim=d_qf-(td-Tj2)*dd_qlimd
    end

    % Compute double S trajectory
    tf=ta+td+tv;
    samples = round((tf-ti)/Ts);
    q=zeros(1,samples);
    d_q=zeros(1,samples);
    dd_q=zeros(1,samples);
    ddd_q=zeros(1,samples);
    ddd_q_min=-ddd_q_max;

    delta_t = tf-ti;
    for i=1:1:samples
        t=tf*(i-1)/(samples-1);
        if t<=Tj1
            q(i)= qi + d_qi*t + ddd_q_max/6*t^3;
            d_q(i)=d_qi + ddd_q_max/2*t^2;
            dd_q(i)=ddd_q_max*t;
            ddd_q(i)=ddd_q_max;
        elseif (t>Tj1) && ( t<ta-Tj1)
            q(i)= qi + d_qi*t + dd_qlima/6*(3*t^2-3*Tj1*t+Tj1^2);
            d_q(i)=d_qi + dd_qlima*(t-Tj1/2);
            dd_q(i)=dd_qlima;
            ddd_q(i)=0;
        elseif (t>ta-Tj1) && (t<=ta )
            q(i)= qi +1.0/2*(d_qlim+d_qi)*ta - d_qlim*(ta-t) - 1.0/6*ddd_q_min*(ta-t)^3;
            d_q(i)=d_qlim + 1.0/2*ddd_q_min*(ta-t)^2;
            dd_q(i)= - ddd_q_min*(ta-t);
            ddd_q(i)=ddd_q_min;
        elseif (t>ta) && (t<=ta+tv )
            q(i)= qi + 1.0/2*(d_qlim+d_qi)*ta + d_qlim*(t-ta);
            d_q(i)=d_qlim;
            dd_q(i)=0;
            ddd_q(i)=0;
        elseif (t>tf-td) && (t<=tf-td+Tj2 )
            q(i)= qf - (d_qlim+d_qf)*td/2 + d_qlim*(t-delta_t+td) - ddd_q_max/6*(t-delta_t+td)^3;
            d_q(i)=d_qlim - ddd_q_max/2*(t-delta_t+td)^2;
            dd_q(i)= - ddd_q_max*(t-delta_t+td);
            ddd_q(i)= - ddd_q_max;
        elseif (t>tf-td+Tj2) && (t<=tf-Tj2 )
            q(i)= qf - (d_qlim+d_qf)*td/2 + d_qlim(t-delta_t-td) + dd_qlimd/6*(3*(t-delta_t+td)^2-3*Tj2*(t-delta_t+td)+Tj2^2);
            d_q(i)= d_qlim + dd_qlimd*(t-delta_t+td-Tj2/2);
            dd_q(i)=-ddd_q_max*Tj2;
            ddd_q(i)=0;
        elseif (t>tf-Tj2) && (t<=tf )
            q(i)= qf - d_qf*(delta_t-t) - ddd_q_max/6*(delta_t-t)^3;
            d_q(i)=d_qf + ddd_q_max/2*(delta_t-t)^2;
            dd_q(i)=- ddd_q_max*(delta_t-t);
            ddd_q(i)=ddd_q_max;
        end
        q(i) = sigma*q(i);
        d_q(i) = sigma*d_q(i);
        dd_q(i) = sigma*dd_q(i);
        ddd_q(i) = sigma*ddd_q(i);
    end
end

