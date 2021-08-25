function [P, DP, DDP, DDDP, T,  p, dp, ddp] = circularMotion(pi, pf, Ts, c, P, DP, DDP,DDDP, T, flag)
    
    

    if size(T)== 0 
        s_tmp  = 0: Ts: norm(pf - pi);
    else
        s_tmp = T(length(T)) + Ts: Ts: T(length(T)) + norm(pf - pi);
    end
    
    if flag == true
        
        userResponse = input('Enter another P. Write the numbers separated by ,\n', 's')
	
        userResponse = strrep(userResponse, ',', ' ');
	
        p_new = sscanf(userResponse, '%f')
        
        [P, DP, DDP, DDDP, T] = circularMotion(pi, p_new, Ts, c, P, DP, DDP, DDDP, T, false);
        [P, DP, DDP, DDDP, T] = circularMotion(p_new, pf, Ts, c, P, DP, DDP, DDDP, T, false);
    
    else
   
        %s = 0: Ts: norm(pf - pi);

        %r = cross(c - pf, c - pi);
        %rho = norm(pi - c);
    
        rho = norm(pi-c);
        r = cross(c-pf, c-pi); 
        arc_angle = vecangle(c-pf, c-pi, r); 
        plength = rho*arc_angle;
    
        s = 0: Ts: plength;
    
        p_prime = [rho*cos((s/rho));rho*sin((s/rho));zeros(1,length(s))];
    
        x_prime = (pi - c)/rho;
        z_prime = r/norm(r);
        y_prime = cross(x_prime, z_prime);
    
        R = [x_prime y_prime z_prime];
    
        p = c + R*p_prime;
        dp = R*[    -sin(s/rho);
                    cos(s/rho);
                    zeros(1,length(s))                       ];
           
        ddp = R*[   -(1/rho)*cos(s/rho);
                    -(1/rho)*sin(s/rho);
                    zeros(1,length(s))                      ];
                
        dddp = R*[  sin(s/rho)/rho^2; 
                    -cos(s/rho)/rho^2; 
                    zeros(1, length(s))                    ]; 
             
        P = [P, p];
        DP = [DP, dp];
        DDP = [DDP, ddp];
        DDDP = [DDDP,dddp];
        T = [T, s_tmp];
    end   
   
end
           
           

