function [relpos, angles, Inn] = get_angles(relpos,rot_angles,m,r,rotation)
%% FUNCTION to randomly rotate an nmer about its com
% INPUTS:
% OUTPUTS:

if rotation
    % make matrices
    Rx = zeros(3,3);
    Ry = zeros(3,3);
    Rz = zeros(3,3);

    % get elements
    a1 = rot_angles(1);
    a2 = rot_angles(2);
    a3 = rot_angles(3);
    c1 = cos(a1);
    c2 = cos(a2);
    c3 = cos(a3);
    s1 = sin(a1);
    s2 = sin(a2);
    s3 = sin(a3);

    % fill matrices

    % Rx
    % first row
    Rx(1,1) = 1; Rx(1,2) = 0; Rx(1,3) = 0;
    % second row
    Rx(2,1) = 0; Rx(2,2) = c1; Rx(2,3) = -s1;
    % third row
    Rx(3,1) = 0; Rx(3,2) = s1; Rx(3,3) = c1;

    % Ry
    % first row
    Ry(1,1) = c2; Ry(1,2) = 0; Ry(1,3) = s2;
    % second row
    Ry(2,1) = 0; Ry(2,2) = 1; Ry(2,3) = 0;
    % third row
    Ry(3,1) = -s2; Ry(3,2) = 0; Ry(3,3) = c2;

    % Rz
    % first row
    Rz(1,1) = c3; Rz(1,2) = -s3; Rz(1,3) = 0;
    % second row
    Rz(2,1) = s3; Rz(2,2) = c3; Rz(2,3) = 0;
    % third row
    Rz(3,1) = 0; Rz(3,2) = 0; Rz(3,3) = 1;

    % euler
    euler = Rx*Ry*Rz;

    NA = size(relpos,1);
    for nn = 1:NA
        rel = relpos(nn,:)';    
        rel = euler*rel;
        relpos(nn,:) = rel';
    end
end

% calculate moment of inertia tensor
Itmp = zeros(3,3);
for i = 1:3
    for j = i:3
        if i == j
            Itmp(i,j) = sum(2/5*m.*r.^2 + m.*(sum(relpos.^2,2)-relpos(:,i).^2));
        else
            Itmp(i,j) = sum(-m.*relpos(:,i).*relpos(:,j));
        end
        Itmp(j,i) = Itmp(i,j);
    end
end

% calculate evals & evecs
[pax, vals_tmp] = eig(Itmp);
Inn = diag(vals_tmp);

% measure euler angles using documentation (see function), output angles
[a1,a2,a3] = meas_euler_angles(pax);
angles = [a1 a2 a3];


end