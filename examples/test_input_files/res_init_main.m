function res_init_main(N,seed,odir,rpdir)

%% FUNCTION to generate random initial configurations of rcp residues
% to be read into sunion_rcp_master.cpp

L = 10*N;
phi0 = 0.1;


%% Get residues
resList = {'Ala','Ile','Leu','Met','Phe','Val'};
NPR = length(resList);

% get residue ratios
fid = fopen([rpdir 'residue_ratios.txt']);
ratiodata = textscan(fid,'%s %f');
ratiolist = ratiodata{1};
ratios = ratiodata{2};
if fid == -1 || isempty(ratios)
    error('ratio file not opened');
end
fclose(fid);

% get atomic radii
fid = fopen([rpdir 'atomic_radii.txt']);
radiidata = textscan(fid,'%s %s %f');
if fid == -1 || isempty(radiidata)
    error('radii file not opened');
end
fclose(fid);

all_radii = cell(NPR,1);
for rn = 1:NPR
    res = resList{rn};
    inds = strcmpi(res,radiidata{1});
    rads = radiidata{3};
    all_radii{rn} = rads(inds);
end

% get fraction of residue types in packing
occupn = round(ratios*N);
occupn_test = occupn;
NT = sum(occupn);
km = 10000;
k = 0;
while NT ~= N && k < km
    k = k + 1;
    occupn_test(randi(NPR)) = occupn(randi(NPR)) + (1-2*randi([0,1]));
    NT = sum(occupn_test(occupn_test>0));
    if NT == N
        occupn = occupn_test;
    end
end

%% Print info for seed

fprintf('Printing info for seed %d\n',seed);
fprintf('Inputting from residue info in %s\n',rpdir);
fprintf('Outputting to %s\n',odir);

% get output file name
outf = [odir 'res_input_N' num2str(N) '_seed' num2str(seed) '.dat'];
outfid = fopen(outf,'w');
format_header = '# id, adiam, ax, ay, az, pdiam, px, py, pz, Ixx, Iyy, Izz, alpha, beta, gamma, M, Na #';

fprintf(outfid,'0\n');
fprintf(outfid,'%d\n',N);
fprintf(outfid,'%f\n',L);
fprintf(outfid,'%f\n',phi0);
fprintf(outfid,'%s\n',format_header);

% id counter
resID = 0;

% monte-carlo initialization
SboxALL = zeros(N,1);
RboxALL = zeros(N,3);

% loop over residues
for rn = 1:NPR        
    rads = all_radii{rn};    
    if occupn == 0
        continue;
    end
    rp = resList{rn};    
    fprintf('* printing residues %s\n', rp);
    fprintf('* rads = %f\n',rads);

    % get occupation number
    rti = strcmpi(ratiolist,rp);
    ni = occupn(rti);

    % get file
    resf = [rpdir rp '_coordinates_whole.txt'];
    fid = fopen(resf);
    rdata = textscan(fid, '%*s %*s %d %f %f %f');
    fclose(fid);

    idALL = rdata{1};
    xALL = rdata{2};
    yALL = rdata{3};
    zALL = rdata{4};
    NALL = idALL(end);   

    % get random inds
    inds = randi(NALL,ni,1);

    % loop over random residues, print info to file
    for r = 1:ni 
        fprintf('** on index %d\n', inds(r));            

        idr = idALL==inds(r);      
        idtmp = r;
        xtmp = xALL(idr);
        ytmp = yALL(idr);
        ztmp = zALL(idr);
        na = length(xtmp);     
        fprintf('na = %d for residue id = %d\n',na,inds(r));

        % get com position
        mi = (4/3)*pi.*rads.^3;
        Msum = sum(mi);
        fprintf('Msum = %f, length(mi) = %d\n',Msum,length(mi));
        xcom = (1/Msum).*sum(mi.*xtmp);
        ycom = (1/Msum).*sum(mi.*ytmp);
        zcom = (1/Msum).*sum(mi.*ztmp);
        R = [xcom ycom zcom];
        rrel = [xtmp ytmp ztmp] - R;                        

        % get initial position in box
        fprintf(' -- getting initial position, mass...');
        Sbox = 2*max(sqrt(sum(rrel(:,1).^2 + rrel(:,2).^2 + rrel(:,3).^2,2))+rads);
        SboxALL(resID+1) = Sbox;

        % give random initial position
        Rbox = L*rand(1,3);
        RboxALL(resID+1,:) = Rbox;
        rbox = rrel + Rbox;
        rbox = mod(rbox,L);                        

        % get residue mass
        dg = 0.25;
        M = get_mass(rrel,rads,Sbox,dg);
        fprintf('.done!\n');

        % get Inn, euler angles
        rot_angles = 2*pi*(2*rand(1,3)-1);
        rotation = 1;
        [rrel, angles, Inn] = get_angles(rrel,rot_angles,mi,rads,rotation);

        % print
        % # id, adiam, ax, ay, az, pdiam, px, py, pz, Ixx, Iyy, Izz, alpha, beta, gamma, M, Na #
        for aa = 1:na
            adiamp = 2*rads(aa);
            xp = rrel(aa,1)+Rbox(1);
            yp = rrel(aa,2)+Rbox(2);
            zp = rrel(aa,3)+Rbox(3);                                
            fprintf(outfid,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n',...
                resID, adiamp, xp, yp, zp, Sbox, Rbox(1), Rbox(2), Rbox(3), ...
                Inn(1), Inn(2), Inn(3), angles(1), angles(2), angles(3), M, na);
        end
        resID = resID + 1;
    end        
end    

end





