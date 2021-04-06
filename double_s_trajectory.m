%% double_s_trajectory.m
%  double_s_trajectory implements double s trajectory taking into account
%  all possible scenario (returning an error in case of invalid
%  constraints). If delta_T is not NaN it will be considered as a special
%  case (fixed duration) and alfa and beta will be considered as well.
%  Output:
%         - returns all trajectories (pos,vel,acc,jerk) and computed tf
function [q,d_q,dd_q,ddd_q,tf] = double_s_trajectory(qi_d,qf_d,d_qi_d,d_qf_d,vd_max,vd_min,ad_max,ad_min,jd_max,jd_min,delta_T,alfa,beta,ti,Ts)
    % Calculate sigma for sign flip
    sigma = sign(qf_d-qi_d);

    % Calculate the correct sign for position and velocity
    qi = sigma * qi_d;
    qf = sigma * qf_d;
    d_qi = sigma * d_qi_d;
    d_qf = sigma * d_qf_d;

    % Calculate the correct values for contraints
    d_qmax = (sigma+1)/2*vd_max + (sigma-1)/2*vd_min;
    % d_qmin = (sigma+1)/2*vd_min + (sigma-1)/2*vd_max;
    dd_qmax = (sigma+1)/2*ad_max + (sigma-1)/2*ad_min;
    % dd_qmin = (sigma+1)/2*ad_min + (sigma-1)/2*ad_max;
    ddd_qmax = (sigma+1)/2*jd_max + (sigma-1)/2*jd_min;
    ddd_qmin = (sigma+1)/2*jd_min + (sigma-1)/2*jd_max;

    % Case 3: Fixed duration provided (with alfa and beta params)
    if not(isnan(delta_T))
      assert(alfa>0 && alfa<=1/2, "alfa params 0<alfa<=1/2")
      assert(beta>0 && beta<=1/2, "beta params 0<beta<=1/2")
      ta = alfa*delta_T;
      td = ta;
      Tj1 = beta*ta;
      Tj2 = Tj1;
      d_qmax = (qf-qi)/((1-alfa)*delta_T);
      dd_qmax = (qf-qi)/(alfa*(1-alfa)*(1-beta)*delta_T^2);
      ddd_qmax = (qf-qi)/(alfa^2*beta*(1-alfa)*(1-beta)*delta_T^3);
      ddd_qmin = -ddd_qmax;
    end
    
    % Feasibility conditions
    cond1 = sqrt(abs(d_qf-d_qi)/ddd_qmax);
    cond2 = dd_qmax/ddd_qmax;
    if cond1 < cond2
        tj = cond1;
    else
        tj = cond2;
    end
    if tj < cond2
        assert((qf-qi)>(tj*(d_qi+d_qf)),"tj < dd_qmax/ddd_qmax: feasibility condition not satisfied");
    elseif tj == cond2
        assert((qf-qi)>1/2*((d_qi+d_qf)*(tj+(d_qf-d_qi)/ddd_qmax)),"tj == dd_qmax/ddd_qmax: feasibility condition not satisfied");
    end

    % We start looking for case 1 since we will be able to move to
    % case 2 in case of tv <= 0
    % Case 1: d_qlim = d_qmax
    if( ((d_qmax - d_qi)*ddd_qmax) < dd_qmax^2 ) 
      % dd_qmax is not reached
      Tj1 = sqrt((d_qmax-d_qi)/ddd_qmax);
      ta = 2*Tj1;
    else
      % dd_qmax is actually reached
      Tj1 = dd_qmax/ddd_qmax;
      ta = Tj1 + (d_qmax - d_qi)/dd_qmax;
    end
    if( ((d_qmax - d_qf)*ddd_qmax) < dd_qmax^2 )
      % dd_qmin is not reached
      Tj2 = sqrt((d_qmax-d_qf)/ddd_qmax);
      td = 2*Tj2;
    else
      % dd_qmin is actually reached
      Tj2 = dd_qmax/ddd_qmax;
      td = Tj2 + (d_qmax-d_qf)/dd_qmax;
    end

    % Velocity time phase
    tv = (qf-qi)/d_qmax - ta/2*(1+d_qi/d_qmax) - td/2*(1+d_qf/d_qmax);

    % Case 2 q_lim < d_qmax
    if(tv < 0)
      % This means that d_qlim < d_qmax
      tv = 0;
      % Tj1=Tj2=Tj
      Tj = dd_qmax/ddd_qmax;
      Tj1 = Tj;
      Tj2 = Tj;
      D = dd_qmax^4/ddd_qmax^2 + 2*(d_qi^2+d_qf^2) + dd_qmax*(4*(qf-qi)-2*dd_qmax/ddd_qmax*(d_qi+d_qf));
      ta = (dd_qmax^2/ddd_qmax - 2*d_qi + sqrt(D))/(2*dd_qmax);
      td = (dd_qmax^2/ddd_qmax - 2*d_qf + sqrt(D))/(2*dd_qmax); 

      assert(ta>=2*Tj, "ta: the maximum (minimum) acceleration is not reached: decrease the value of dd_qmax");
      assert(td>=2*Tj, "td: the maximum (minimum) acceleration is not reached: decrease the value of dd_qmax");
      % if ta<0 or td<0 only one of the acceleration or deceleration phase is necessary
    end
    
    % Just and extra check that comes from nonsense constraints
    assert(isreal(ta), "ta: is not a real number: change limit values");
    assert(isreal(td), "td: is not a real number: change limit values");

    % Computation section
    % Compute acceleration limits
    dd_qlima = ddd_qmax*Tj1;
    dd_qlimd = -ddd_qmax*Tj2;
    d_qlim = d_qi + (ta-Tj1)*dd_qlima; % Same as d_qlim=d_qf-(td-Tj2)*dd_qlimd
        
    % Compute total duration (tf)
    duration=ta+tv+td;
    % Compute final time of the trajectory
    tf=duration+ti;

    % Compute time instances to calculate the trajectory
    time=ti:Ts:tf;
    delta_t = tf-ti;

    % Preallocate just to avoid MATLAB warning
    samples = round((tf-ti)/Ts);
    q_hat=zeros(1,samples); 
    d_qhat=zeros(1,samples);
    dd_qhat=zeros(1,samples); 
    ddd_qhat=zeros(1,samples);
    
    % Preallocate also for the final result
    q=zeros(1,samples); 
    d_q=zeros(1,samples);
    dd_q=zeros(1,samples); 
    ddd_q=zeros(1,samples);
    i = 1;

    % Acually compute the trajectory by splitting each segment and
    % calculate position, velocity, acceleration and jerk
    for t = time
      if( t <= Tj1+ti )
        q_hat(i) = qi + d_qi*(t-ti) + ddd_qmax * (t-ti)^3/6;
        d_qhat(i) = d_qi + ddd_qmax*(t-ti)^2/2;
        dd_qhat(i) = ddd_qmax * (t-ti);
        ddd_qhat(i) = ddd_qmax;
      elseif( (t > Tj1+ti) && (t <= (ta-Tj1)+ti) )
        q_hat(i) = qi + d_qi*(t-ti) + dd_qlima/6*(3*(t-ti)^2-3*Tj1*(t-ti)+Tj1^2);
        d_qhat(i) = d_qi + dd_qlima*((t-ti)-Tj1/2);
        dd_qhat(i) = ddd_qmax*Tj1;
        ddd_qhat(i) = 0;
      elseif ( (t > (ta-Tj1)+ti) && (t <= ta+ti) )
        q_hat(i) = qi + (d_qlim + d_qi)*ta/2 - d_qlim*(ta-(t-ti)) - ddd_qmin*(ta-(t-ti))^3/6;
        d_qhat(i) = d_qlim + ddd_qmin*(ta-(t-ti))^2/2;
        dd_qhat(i) = -ddd_qmin*(ta-(t-ti));
        ddd_qhat(i) = ddd_qmin;
      elseif( (t > ta+ti) && (t <= (ta+tv)+ti) )
        q_hat(i) = qi + (d_qlim+d_qi)*ta/2 + d_qlim*((t-ti)-ta);
        d_qhat(i) = d_qlim;
        dd_qhat(i) = 0;
        ddd_qhat(i) = 0;
      elseif( (t > (ta+tv+ti)) && (t <= (ta+tv+Tj2+ti)) )
        q_hat(i) = qf - (d_qlim+d_qf)*td/2 + d_qlim*((t-ti)-delta_t+td) - ddd_qmax*((t-ti)-delta_t+td)^3/6;
        d_qhat(i) = d_qlim - ddd_qmax*((t-ti)-delta_t+td)^2/2;
        dd_qhat(i) = -ddd_qmax*((t-ti)-delta_t+td);
        ddd_qhat(i) = ddd_qmin;
      elseif( (t > (ta+tv+Tj2+ti)) && (t <= (ta+tv+(td-Tj2+ti))) )
        q_hat(i) = qf-(d_qlim+d_qf)*td/2+d_qlim*((t-ti)-delta_t+td)+dd_qlimd/6*(3*((t-ti)-delta_t+td)^2-3*Tj2*((t-ti)-delta_t+td)+Tj2^2); 
        d_qhat(i) = d_qlim + dd_qlimd*((t-ti)-delta_t+td-Tj2/2);
        dd_qhat(i) = -ddd_qmax*Tj2;
        ddd_qhat(i) = 0;
      elseif( (t > (ta+tv+(td-Tj2+ti))) && (t <= tf) )
        q_hat(i) = qf - d_qf*(delta_t-(t-ti)) - ddd_qmax*(delta_t-(t-ti))^3/6;
        d_qhat(i) = d_qf + ddd_qmax*(delta_t-(t-ti))^2/2;
        dd_qhat(i) = -ddd_qmax*(delta_t-(t-ti));
        ddd_qhat(i) = ddd_qmax;
      end

      q(i) = sigma*q_hat(i);
      d_q(i) = sigma*d_qhat(i);
      dd_q(i) = sigma*dd_qhat(i);
      ddd_q(i) = sigma*ddd_qhat(i);
      i=i+1;
    end
    
    % Print value computed for debug purposes
    fprintf(' ti=%f\n tj1=%f\n ta=%f\n tv=%f\n tj2=%f\n td=%f\n tf=%f\n dd_qlima=%f\n dd_qlimd=%f\n d_qlim=%f\n', ti,Tj1,ta,tv,Tj2,td,tf, sigma*dd_qlima, sigma*dd_qlimd,sigma*d_qlim);
end