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
ddd_qmax=3;
ddd_qmin=ddd_qmax;
ti = 0;
Ts = 0.01;

[q,d_q,dd_q,ddd_q, tf]=double_s_trajectory(qi,qf,d_qi,d_qf,d_qmin,d_qmax,dd_qmin,dd_qmax,ddd_qmin,ddd_qmax,ti,Ts);
t = linspace(ti, tf, (tf-ti)/Ts);
plot_trajectory('Double S Trajectory', t, [q; d_q; dd_q; ddd_q])

function  [q,d_q,dd_q,ddd_q, tf]=double_s_trajectory(qi,qf,d_qi,d_qf,d_qmin,d_qmax,dd_qmin,dd_qmax,ddd_qmin,ddd_qmax, ti,Ts)
    % Case 1: d_qlim = d_qmax
    if (d_qmax-d_qi)*ddd_qmax<dd_qmax^2    
        % dd_qmax is not reached
        Tj1=sqrt((d_qmax-d_qi)/ddd_qmax);
        ta=2*Tj1;
    else                        
        % dd_qmax is actually reached
        Tj1=dd_qmax/ddd_qmax;
        ta=Tj1+(d_qmax-d_qi)/dd_qmax;
    end
    if (d_qmax-d_qf)*ddd_qmax<dd_qmax^2    
        % dd_qmin is not reached
        Tj2=sqrt((d_qmax-d_qf)/ddd_qmax);
        td=2*Tj2;
    else                       
        % dd_qmin is actually reached
        Tj2=dd_qmax/ddd_qmax;
        td=Tj2+(d_qmax-d_qf)/dd_qmax;
    end

    tv=(qf-qi)/d_qmax-ta/2*(1+d_qi/d_qmax)-td/2*(1+d_qf/d_qmax);

    if tv>0
        % The maximum velocity is actually reached
        d_qlim=d_qmax;
        dd_qlima=ddd_qmax*Tj1;
        dd_qlimd=-ddd_qmax*Tj2;
    end

    % Case 2 d_qlim < d_qmax
    if tv<=0 
        % This means that d_qlim < d_qmax
        tv=0;
        Tj1=dd_qmax/ddd_qmax;
        Tj2=Tj1;
        tj=Tj1;

        D=dd_qmax^4/(ddd_qmax^2)+2*(d_qi^2+d_qf^2)+dd_qmax*(4*(qf-qi)-2*(dd_qmax/ddd_qmax)*(d_qi+d_qf));
        ta=(dd_qmax^2/ddd_qmax-2*d_qi+sqrt(D))/(2*dd_qmax);
        td=(dd_qmax^2/ddd_qmax-2*d_qf+sqrt(D))/(2*dd_qmax);
        assert(isreal(ta), "ta is not a real number");
        assert(isreal(td), "td is not a real number");

        assert(ta<2*tj || td<2*tj, "the maximum (minimum) acceleration is not reached: decrease the value of dd_qmax");
        % assert(ta<0 || td<0, "only one of the acceleration or deceleration phase is necessary");
        if (ta<0 || td<0)
            display("only one of the acceleration or deceleration pahse is necessary")
        end
    end
    
    dd_qlima=ddd_qmax*Tj1;
    dd_qlimd=-ddd_qmax*Tj2;
    d_qlim=d_qi+(ta-Tj1)*dd_qlima; % Same as d_qlim=d_qf-(td-Tj2)*dd_qlimd

    % Compute double S trajectory
    sigma = sign(qf-qi);
    qi = sigma * qi;
    qf = sigma * qf;
    d_qi = sigma * d_qi;
    d_qf = sigma * d_qf;

    d_qmax = (sigma+1)/2*d_qmax + (sigma-1)/2*d_qmin;
%     d_qmin = (sigma+1)/2*d_qmin + (sigma-1)/2*d_qmax;
    dd_qmax = (sigma+1)/2*dd_qmax + (sigma-1)/2*dd_qmin;
%     dd_qmin = (sigma+1)/2*dd_qmin + (sigma-1)/2*dd_qmax;
    ddd_qmax = (sigma+1)/2*ddd_qmax + (sigma-1)/2*ddd_qmin;
%     ddd_qmin = (sigma+1)/2*ddd_qmin + (sigma-1)/2*ddd_qmax;

    tf=ta+td+tv+ti;
    samples = round((tf-ti)/Ts);
    q=zeros(1,samples);
    d_q=zeros(1,samples);
    dd_q=zeros(1,samples);
    ddd_q=zeros(1,samples);
    ddd_qmin=-ddd_qmax;

    delta_t = tf-ti;
    for i=1:1:samples
        t=tf*(i-1)/(samples-1);
        if t<=Tj1+ti
            q(i)= qi + d_qi*(t-ti) + ddd_qmax/6*(t-ti)^3;
            d_q(i)=d_qi + ddd_qmax/2*(t-ti)^2;
            dd_q(i)=ddd_qmax*(t-ti);
            ddd_q(i)=ddd_qmax;
        elseif (t>Tj1+ti) && (t<ta-Tj1)
            q(i)= qi + d_qi*(t-ti) + dd_qlima/6*(3*(t-ti)^2-3*Tj1*(t-ti)+Tj1^2);
            d_q(i)=d_qi + dd_qlima*((t-ti)-Tj1/2);
            dd_q(i)=dd_qlima;
            ddd_q(i)=0;
        elseif (t>ta-Tj1) && (t<=ta+ti)
            q(i)= qi +1.0/2*(d_qlim+d_qi)*ta - d_qlim*(ta-(t-ti)) - 1.0/6*ddd_qmin*(ta-(t-ti))^3;
            d_q(i)=d_qlim + 1.0/2*ddd_qmin*(ta-(t-ti))^2;
            dd_q(i)= - ddd_qmin*(ta-(t-ti));
            ddd_q(i)=ddd_qmin;
        elseif (t>ta+ti) && (t<=ta+tv+ti)
            q(i)= qi + 1.0/2*(d_qlim+d_qi)*ta + d_qlim*((t-ti)-ta);
            d_q(i)=d_qlim;
            dd_q(i)=0;
            ddd_q(i)=0;
        elseif (t>tf-td) && (t<=tf-td+Tj2 )
            q(i)= qf - (d_qlim+d_qf)*td/2 + d_qlim*((t-ti)-delta_t+td) - ddd_qmax/6*((t-ti)-delta_t+td)^3;
            d_q(i)=d_qlim - ddd_qmax/2*((t-ti)-delta_t+td)^2;
            dd_q(i)= - ddd_qmax*((t-ti)-delta_t+td);
            ddd_q(i)= - ddd_qmax;
        elseif (t>tf-td+Tj2) && (t<=tf-Tj2 )
            q(i)= qf - (d_qlim+d_qf)*td/2 + d_qlim((t-ti)-delta_t-td) + dd_qlimd/6*(3*((t-ti)-delta_t+td)^2-3*Tj2*((t-ti)-delta_t+td)+Tj2^2);
            d_q(i)= d_qlim + dd_qlimd*((t-ti)-delta_t+td-Tj2/2);
            dd_q(i)=-ddd_qmax*Tj2;
            ddd_q(i)=0;
        elseif (t>tf-Tj2) && (t-ti<=tf)
            q(i)= qf - d_qf*(delta_t-(t-ti)) - ddd_qmax/6*(delta_t-(t-ti))^3;
            d_q(i)=d_qf + ddd_qmax/2*(delta_t-(t-ti))^2;
            dd_q(i)=- ddd_qmax*(delta_t-(t-ti));
            ddd_q(i)=ddd_qmax;
        end
        q(i) = sigma*q(i);
        d_q(i) = sigma*d_q(i);
        dd_q(i) = sigma*dd_q(i);
        ddd_q(i) = sigma*ddd_q(i);
    end
end

