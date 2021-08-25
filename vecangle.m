function angle = vecangle(v1, v2, normal)

    toRad = 2*pi/360; 
    xprod = cross(v1,v2);
    c = sign(dot(xprod,normal)) * norm(xprod); 
    angle = atan2d(c,dot(v1,v2)) * toRad; 
    
end 

