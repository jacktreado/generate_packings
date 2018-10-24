function [alpha,beta,gamma] = meas_euler_angles(pax)
% Based on euler angle document in 
% ~/Jamming/literature/proteins/molecular_dynamics/

% relabel Principal Axes P1,P2,P3
P1 = pax(:,1);
P2 = pax(:,2);
P3 = pax(:,3);

P1x = P1(1);
P1y = P1(2);
P1z = P1(3);

P2x = P2(1);
P2y = P2(2);
P2z = P2(3);

P3x = P3(1);
P3y = P3(2);
P3z = P3(3);

if (abs(sum(P3.*[0; 0; 1])-1) < 1e-16)
    alpha = 0;
    if (P3z > 0)
        beta = 0;
    else
        beta = pi;
    end
    gamma = -atan2(P1y,P1x);
else
    alpha = atan2(P3x,-P3y);
    beta = acos(P3z);
    gamma = atan2(P1z,P2z);
end

end