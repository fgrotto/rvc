function [P, DP, DDP, DDDP, T] = pipeline(Pg, C, tag, motion, Ts)

    P = [];
    DP = [];
    DDP = [];
    DDDP = [];
    T = [];
    j = 1;
    
    for i= 1:length(Pg)-1
        
        if motion(i) == 1
            
            [P, DP, DDP, DDDP, T] = linearMotion(Pg(:,i), Pg(:,i+1), Ts, P, DP, DDP,DDDP, T);
        
        else
            
            if tag(j) == 0
                [P, DP, DDP, DDDP, T] = circularMotion(Pg(:,i), Pg(:,i+1), Ts, C(:,j), P, DP, DDP, DDDP, T, true);
            else
                [P, DP, DDP, DDDP, T] = circularMotion(Pg(:,i), Pg(:,i+1), Ts, C(:,j), P, DP, DDP, DDDP, T, false);
            end
        
            j = j + 1;
        end
     
    end

            
end

