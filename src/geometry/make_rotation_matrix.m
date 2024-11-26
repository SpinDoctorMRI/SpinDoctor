function R = make_rotation_matrix(direction)
    z = [0,0,1];
    direction = direction/norm(direction);
    if dot(z, direction) == -1
        R = [1,0,0;0,-1,0 ;0,0,-1];
    elseif dot(z, direction) == 1
        R = eye(3);
    else
        v = cross(z,direction);
        s = norm(v);
        c  = dot(z,direction);
        v_plus = [0, - v(3) ,  v(2);v(3), 0 , -v(1);-v(2), v(1), 0 ];
        if s > 0.15
            R = eye(3) + v_plus + ((1 -c)/s^2)*v_plus^2;
        else
            R = eye(3) + v_plus + (1/2 - s^2/24)*v_plus^2;
        end
    end

end