%% multipoint_trajectory.m
%  It calculates a multipoint trajectory given the desired pairs (qk, tk)
%  It used the given d_qi and d_qf provided
%  It returns:
%         - The position, velocitiy and acceleration profiles
function [q,d_q,dd_q] = multipoint_trajectory(q_d, t_d, d_qi, d_qf, Ts)
    assert(length(q_d) == length(t_d), "they should be pairs of (q_k,t_k)");
    n = length(q_d);

    % Calculate v vector according to the Euler approximation
    for k=2:n
      v(k) = (q_d(k)-q_d(k-1))/(t_d(k)-t_d(k-1));   
    end

    % Check that the time is sorted correctly @TODO

    % Calculate velocity in each path point pay attention to indexes
    % since in matlab the first position is 1!
    v(1) = d_qi;
    for k=2:n-1
      if (sign(v(k)) == sign(v(k+1)))
        v(k) = (v(k)+v(k+1))/2;
      else
        v(k) = 0;
      end
    end
    v(n) = d_qf;

    i = 1;
    % Calculate each trajectory segment pi(k) and sum the position,
    % velocity and acceleration profile obtained
    for k=1:n-1
      Tk = t_d(k+1)-t_d(k);
      a0 = q_d(k);
      a1 = v(k);
      a2 = 1/Tk*(3*(q_d(k+1)-q_d(k))/Tk - 2*v(k) - v(k+1));
      a3 = 1/Tk^2*(2*(q_d(k)-q_d(k+1))/Tk+v(k)+v(k+1));

      % Pay attention to the time and avoid sample repetitions
      time=t_d(k):Ts:t_d(k+1);
      if k == n-1
         time = time(1:end);
      else
        time = time(1:end-1);
      end
      
      % Add the calculated sample to the final position, velocities and accelerations
      for t = time
        q(i) = a3*(t-t_d(k))^3+a2*(t-t_d(k))^2+a1*(t-t_d(k))+a0;
        d_q(i) = 3*a3*(t-t_d(k))^2+2*a2*(t-t_d(k))+a1;
        dd_q(i) = 6*a3*(t-t_d(k))+2*a2;
        i=i+1;
      end
    end
end
