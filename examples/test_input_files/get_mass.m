function m = get_mass(relpos,r,s,g)
%% Function to calculate mass of given residue

% put in center of bounding box
N = size(relpos,1);
L = s + 2*g;
rcom_c = [L/2 L/2 L/2];
pos_c = relpos + rcom_c;

% monte carlo integration
mv = 0;

% make grid
G = round(L/g);
g = L/G;
[gx, gy, gz] = ndgrid(0:g:L-g);


for nx = 1:G
    for ny = 1:G
        for nz = 1:G
            ptx = gx(nx,ny,nz);
            pty = gy(nx,ny,nz);
            ptz = gz(nx,ny,nz);
            pt = [ptx, pty, ptz];
            outside = 1;
            for n = 1:N
                apos = pos_c(n,:);
                rad = r(n);
                d = apos-pt;
                dist = sqrt(sum(d.^2));
                if dist < rad
                    outside = 0;
                    break;                    
                end
            end
            
            % if not in residue, get outside vol
            if outside
                mv = mv + g^3;
            end
        end
    end
end

m = L^3-mv;