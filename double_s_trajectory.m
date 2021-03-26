close all; clear all;

qi_d = 0;
qf_d = 2;
d_qi_d = 0;
d_qf_d = 0;
ti = 5;
Ts = 0.001;

vd_max = 4;  vd_min = -vd_max;
ad_max = 20;  ad_min = -ad_max;
jd_max = 60;  jd_min = -jd_max;

% Calculate sigma for sign flip
sigma = sign(qf_d-qi_d);

% Calculate the correct sign for position and velocity
qi = sigma * qi_d;
qf = sigma * qf_d;
d_qi = sigma * d_qi_d;
d_qf = sigma * d_qf_d;

% Calculate the correct values for contrains
d_qmax = (sigma+1)/2*vd_max + (sigma-1)/2*vd_min;
d_qmin = (sigma+1)/2*vd_min + (sigma-1)/2*vd_max;
dd_qmax = (sigma+1)/2*ad_max + (sigma-1)/2*ad_min;
dd_qmin = (sigma+1)/2*ad_min + (sigma-1)/2*ad_max;
ddd_qmax = (sigma+1)/2*jd_max + (sigma-1)/2*jd_min;
ddd_qmin = (sigma+1)/2*jd_min + (sigma-1)/2*jd_max;

% We start looking for case 1 since we will be able to move to
% case 2 in case of tv <= 0
% Case 1: d_qlim = d_qmax
if( ((d_qmax - d_qi)*ddd_qmax) < dd_qmax^2 ) 
  % dd_qmax is not reached
  Tj1 = sqrt((d_qmax-d_qi)/ddd_qmax);
  Ta = 2*Tj1;
else
  % dd_qmax is actually reached
  Tj1 = dd_qmax/ddd_qmax;
  Ta = Tj1 + (d_qmax - d_qi)/dd_qmax;
end
if( ((d_qmax - d_qf)*ddd_qmax) < dd_qmax^2 )
  % dd_qmin is not reached
  Tj2 = sqrt((d_qmax-d_qf)/ddd_qmax);
  Td = 2*Tj2;
else
  % dd_qmin is actually reached
  Tj2 = dd_qmax/ddd_qmax;
  Td = Tj2 + (d_qmax-d_qf)/dd_qmax;
end

% Velocity time phase
tv = (qf-qi)/d_qmax - Ta/2*(1+d_qi/d_qmax) - Td/2*(1+d_qf/d_qmax);
if( tv < 0 )
  tv = 0; % must be set to 0
    
  % Reset the loop counter
  count = 0;
  
  % Do this loop until the break condition holds
  for gamma = 1:-0.001:0
    
    dd_qmax = gamma * dd_qmax;
    dd_qmin = gamma * dd_qmin;

    % d_qmax is not reached => check other things
    Tj1 = dd_qmax/ddd_qmax;
    Tj2 = dd_qmax/ddd_qmax;
    Tj = dd_qmax/ddd_qmax;
    delta = dd_qmax^4/ddd_qmax^2 + 2*(d_qi^2+d_qf^2) + dd_qmax*(4*(qf-qi)-2*dd_qmax/ddd_qmax*(d_qi+d_qf));
    Ta = (dd_qmax^2/ddd_qmax - 2*d_qi + sqrt(delta))/(2*dd_qmax);
    Td = (dd_qmax^2/ddd_qmax - 2*d_qf + sqrt(delta))/(2*dd_qmax);    
    
    if( Ta < 0 )   
      Ta = 0;
      Td = 2*(qf-qi)/(d_qf+d_qi);
      Tj2 = (ddd_qmax*(qf-qi) - sqrt(ddd_qmax*(ddd_qmax*(qf-qi)^2+(d_qf+d_qi)^2*(d_qf-d_qi))))/(ddd_qmax*(d_qf+d_qi));
    elseif( Td < 0 )
      Td = 0;
      Ta = 2*(qf-qi)/(d_qf+d_qi);
      Tj1 = (ddd_qmax*(qf-qi) - sqrt(ddd_qmax*(ddd_qmax*(qf-qi)^2-(d_qf+d_qi)^2*(d_qf-d_qi))))/(ddd_qmax*(d_qf+d_qi));
    else
      if( (Ta > 2*Tj) && (Td > 2*Tj) )
        break;
      else
        count = count + 1;
      end
    end
  end
end

dd_qlima = ddd_qmax * Tj1;
dd_qlimd = -ddd_qmax*Tj2;
d_qlim = d_qi + (Ta-Tj1)*dd_qlima;

tf=Ta+tv+Td;
tf=tf+ti;

time=ti:Ts:tf;
delta_t = tf-ti;

% Preallocate just to avoid MATLAB warning
samples = round((tf-ti)/Ts);
q=zeros(1,samples);
qp=zeros(1,samples);
qpp=zeros(1,samples);
qppp=zeros(1,samples);
i = 1;

for t = time
  if( t <= Tj1+ti )
    q(i) = qi + d_qi*(t-ti) + ddd_qmax * (t-ti)^3/6;
    qp(i) = d_qi + ddd_qmax*(t-ti)^2/2;
    qpp(i) = ddd_qmax * (t-ti);
    qppp(i) = ddd_qmax;
  elseif( (t > Tj1+ti) && (t <= (Ta-Tj1)+ti) )
    q(i) = qi + d_qi*(t-ti) + dd_qlima/6*(3*(t-ti)^2-3*Tj1*(t-ti)+Tj1^2);
    qp(i) = d_qi + dd_qlima*((t-ti)-Tj1/2);
    qpp(i) = ddd_qmax*Tj1;
    qppp(i) = 0;
  elseif ( (t > (Ta-Tj1)+ti) && (t <= Ta+ti) )
    q(i) = qi + (d_qlim + d_qi)*Ta/2 - d_qlim*(Ta-(t-ti)) - ddd_qmin*(Ta-(t-ti))^3/6;
    qp(i) = d_qlim + ddd_qmin*(Ta-(t-ti))^2/2;
    qpp(i) = -ddd_qmin*(Ta-(t-ti));
    qppp(i) = ddd_qmin;
  elseif( (t > Ta+ti) && (t <= (Ta+tv)+ti) )
    q(i) = qi + (d_qlim+d_qi)*Ta/2 + d_qlim*((t-ti)-Ta);
    qp(i) = d_qlim;
    qpp(i) = 0;
    qppp(i) = 0;
  elseif( (t > (Ta+tv+ti)) && (t <= (Ta+tv+Tj1+ti)) )
    q(i) = qf - (d_qlim+d_qf)*Td/2 + d_qlim*((t-ti)-delta_t+Td) - ddd_qmax*((t-ti)-delta_t+Td)^3/6;
    qp(i) = d_qlim - ddd_qmax*((t-ti)-delta_t+Td)^2/2;
    qpp(i) = -ddd_qmax*((t-ti)-delta_t+Td);
    qppp(i) = ddd_qmin;
  elseif( (t > (Ta+tv+Tj2+ti)) && (t <= (Ta+tv+(Td-Tj2+ti))) )
    q(i) = qf-(d_qlim+d_qf)*Td/2+d_qlim*((t-ti)-delta_t+Td)+dd_qlimd/6*(3*((t-ti)-delta_t+Td)^2-3*Tj2*((t-ti)-delta_t+Td)+Tj2^2); 
    qp(i) = d_qlim + dd_qlimd*((t-ti)-delta_t+Td-Tj2/2);
    qpp(i) = -ddd_qmax*Tj2;
    qppp(i) = 0;
  elseif( (t > (Ta+tv+(Td-Tj2+ti))) && (t <= tf) )
    q(i) = qf - d_qf*(delta_t-(t-ti)) - ddd_qmax*(delta_t-(t-ti))^3/6;
    qp(i) = d_qf + ddd_qmax*(delta_t-(t-ti))^2/2;
    qpp(i) = -ddd_qmax*(delta_t-(t-ti));
    qppp(i) = ddd_qmax;
  end
  
  qd(i) = sigma*q(i);
  qdp(i) = sigma*qp(i);
  qdpp(i) = sigma*qpp(i);
  qdppp(i) = sigma*qppp(i);
  i=i+1;
end

t = linspace(ti, tf, (tf-ti)/Ts);
t(end+1)=tf;
plot_trajectory('Double S Trajectory', t, [qd; qdp; qdpp; qdppp])