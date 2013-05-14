%% Monte Carlo simulation
%for measuring photons detected with sphere over tissue
%used to calculate scattering losses as a function of optical properties

% Last editted May 13, 2013 at 13:00 by Diana Glennie.

%% MAIN PROGRAM
function [] = mc_tissue_sphere_mega
global g g2 nrel_at nrel_ta reflco rusnum rusfrac kftn israd portrad midz coshalf surfrad mua musp mut albedo stotal stotalbin

% seed the random number generator based on the current time
rng('shuffle');

disp(['The current simulation began at',num2str(clock)])
% read XLSX file containing list of MC sims to run & related data
[file, nruns] = read_list_sims;

% for each simulation,
for currun = 1:nruns %first run starts at 1 (Sheet1)
            %disp(['Currently processing run #',num2str(currun)])
    [g g2 nrel_at nrel_ta reflco rusnum rusfrac kftn israd portrad midz coshalf surfrad mua musp mut albedo] = read_param(currun, file); %read in simulation parameters
    
    stotalbin = zeros(10001,1);
    
    
    % for each group of one thousand photons,
    for kftncount = 1:kftn
        h = waitbar(0,['Currently working on Run ',num2str(currun),'/', num2str(nruns),', Group ', num2str(kftncount),'/',num2str(kftn)],'Name',file);
        % for each photon,
        for ftncount = 1:1000
            waitbar(ftncount/1000)
%              disp(ftncount)
            [pos, dir, nrus, ftnwt, tissuesphere, stotal] = ftnini; % initialize photon
            while ftnwt ~= 0
                if tissuesphere == 1
                    [pos, dir, ftnwt, tissuesphere] = sphere(pos, dir, ftnwt, tissuesphere);
                elseif tissuesphere == -1
                    [pos, dir, ftnwt, tissuesphere, nrus] = tissue(pos, dir, ftnwt, tissuesphere, nrus);
                else
                    disp('Something has gone terribly wrong (tissuesphere =/ +/-1).')
                    ftnwt = 0;
                end

               if ftnwt < 1e-16 % ftnwt too small
                    ftnwt = 0;
                    %disp('Too small')
               end 
                
                if (pos(3) < 0) && (sqrt(pos(1)^2 + pos(2)^2 + pos(3)^2) >= 100) % too far away
                    ftnwt = 0;
                    %radius = sqrt(pos(1)^2 + pos(2)^2 + pos(3)^2)
                    %disp('Too far away')
                end
                %if nrus >= rusnum % too many scatters
                %    [ftnwt,nrus] = rusrul(ftnwt,nrus);
                %end
            end
        end
        
        str = num2str(currun);
        runloc = strcat('Sheet',str);
        xlswrite(file,stotalbin,runloc,'D1')% print detectwt to current run's sheet
        
        close(h)
    end
    
end
disp(['The current run was completed at',num2str(clock)])
end

%% Read in location of data and number of simulations to run
function [file, nruns] = read_list_sims % VERIFIED

[filename, pathname] = uigetfile('*.xlsx', 'Select the XLSX file containing the MC data to run');
file = strcat(pathname,filename);
nruns = xlsread(file,'Sheet0');

end

%% Read in parameters for the current simulation
function [g g2 nrel_at nrel_ta reflco rusnum rusfrac kftn israd portrad midz coshalf surfrad mua musp mut albedo] = read_param(currun, file) % VERIFIED

    str = num2str(currun);
    runloc = strcat('Sheet',str);
    [num_data, text_data] = xlsread(file,runloc);
    
    g = num_data(1); % anisotropy parameter
    g2 = 1- g^2;
    nrel_at = num_data(2); % relative index of refraction (air to tissue)
    nrel_ta = 1/nrel_at;
    reflco = num_data(3); % reflection coefficient of sphere
    rusnum = num_data(4); % # of scatters before Russian roulette
    rusfrac = num_data(5); % survival fraction applied to successful roulettes
    kftn = num_data(6); % number of groups of thousands of photons
    israd = num_data(7); % radius of integrating sphere (mm)
    portrad = num_data(8); % radius of detection port (mm)
    midz = sqrt(israd^2-portrad^2);
    coshalf = num_data(9); % cosine of the fiber half angle
    surfrad = (sqrt(1-coshalf^2)/coshalf)*(israd + midz);
    mua = num_data(10); % absorption coefficient (mm^-1)
    musp = num_data(11); % reduced scatter coefficient (mm^-1)
    mut = mua + musp/(1-g);
    albedo = (mut-mua)/mut; % transport/scatter albedo
    
end

%% Initialize photon
function [pos, dir, nrus, ftnwt, tissuesphere, stotal] = ftnini % VERIFIED
global israd midz coshalf


pos = [israd 0 midz];

costheta = coshalf + (1 - coshalf)*rand;
sintheta = sqrt(1 - costheta^2);

cp = 2;
while ((cp >= 1) || (cp == 0));
   p = 2*rand - 1;
   q = 2*rand - 1;
   cp = q^2 + p^2;
end
cp = sqrt(cp);
sinphi = q/cp;
cosphi = p/cp;

dir = [-costheta sintheta*sinphi sintheta*cosphi];

%ntissuescat = 0;% number of scatters in tissue for scoring purposes
%nspherescat = 0;% number of scatters in sphere for scoring purposes
nrus = 0; % number of scatters for Russian roulette
ftnwt = 1;
tissuesphere = 1; % if =1 then in sphere, if equal -1 then in tissue
stotal = 0;

end

%% Russian roulette
function [ftnwt,nrus] = rusrul(ftnwt,nrus)
global rusfract

    if (rand <= rusfract)
        nrus = -10*nrus;
        ftnwt = ftnwt/rusfrac;
    else
        ftnwt = 0;
        disp('Roulette killed.')
    end
    
end

%% Map how photon moves while in sphere
function [pos, dir, ftnwt, tissuesphere] = sphere(pos, dir, ftnwt, tissuesphere)
global nrel_at nrel_ta reflco israd midz surfrad stotalbin stotal

    oldpos = pos;
    pos = raytrace(pos, dir); % raytrace to new position on sphere
    
    if pos(3) <= 0 % then have gone into tissue
        % backtrack to tissue surface
        pathlength = -pos(3)/dir(3);
        pos = adjust_pos(pos, dir, pathlength);
        
        % at tissue surface, check to see if photon reflects or transmits:
        bscwt = fresnel(nrel_at, dir);
        
        if rand <= bscwt % reflects off surface
            dir(3) = -dir(3);
        else % successfully transmits into tissue
            dir = refract(nrel_ta, dir);
            tissuesphere = -1;
        end
        
    else % photon is still in sphere
        if (sqrt(oldpos(1)^2+oldpos(2)^2) < surfrad) && (abs(oldpos(3)) < 1) && (pos(3) > (midz + sqrt(israd^2 - 2.5^2)))
%  disp('Photon scored')            
            if stotal <= 0
                stotalbin(1) = stotalbin(1) + ftnwt;
            elseif stotal > 1000
                stotalbin(10001) = stotalbin(10001) + ftnwt;
            else
                bin = ceil(stotal*10 + 1); % find bin
                stotalbin(bin) = stotalbin(bin) + ftnwt; % add total pathlength of photon to bin
            end
            ftnwt = 0;
            
        else % is reflected off sphere wall
            ftnwt = ftnwt*reflco; % adjust photon weight
            dir = sphere_scat(pos); % sample for new direction vector
        end
    end

end

%% Calculates the pathlength and then position along the current direction vector to the next spot on the sphere wall
function newpos = raytrace(pos, dir) % VERIFIED
global israd midz
    

    a = 1;
    b = 2*(pos(1)*dir(1) + pos(2)*dir(2) + (pos(3) - midz)*dir(3));
    c = pos(1)^2 + pos(2)^2 + (pos(3) - midz)^2 - israd^2;
    
    if (b^2-4*a*c) < 0
        pathlength = (-b)/(2*a);
    else
        pathlength = (- b + sqrt(b^2 - 4*a*c))/(2*a);
    end
    
    newpos = adjust_pos(pos, dir, pathlength);
end

%% Move from old position along direction vector for determined length to new position
function newpos = adjust_pos(oldpos, dir, pathlength)
    newpos(1) = oldpos(1) + pathlength*dir(1);
    newpos(2) = oldpos(2) + pathlength*dir(2);
    newpos(3) = oldpos(3) + pathlength*dir(3);

end

%% Calculate the Fresnel reflection and transmission coefficients
function bscwt = fresnel(nrel, dir)
    if nrel == 1
        bscwt = 0;
    else
        nrsq = nrel^2;
        snisq = 1 - dir(3)^2;
        snrsq = snisq/nrsq;
        if snrsq >= 1
            bscwt = 1;
        else
            ns = sqrt(nrsq - snisq);
            np = abs(nrsq*dir(3));
            rs = (abs(dir(3)) - ns)/(abs(dir(3)) + ns);
            rp = (ns - np)/(ns + np);
            bscwt = (rs^2 + rp^2)/2;
        end
    end

end

%% Changes direction vector based on angle of refraction
function dir = refract(nrel, dir)
    dir(1) = nrel*dir(1);
    dir(2) = nrel*dir(2);
    dir(3) = sqrt(1-nrel^2*(1-dir(3)^2))*(dir(3)/abs(dir(3)));

end

%% Scatter inside sphere
function dir = sphere_scat(pos)
global midz

    % point toward center of sphere
    pathlength = sqrt(pos(1)^2 + pos(2)^2 + (midz-pos(3))^2);
    dir(1) = -pos(1)/pathlength;
    dir(2) = -pos(2)/pathlength;
    dir(3) = (midz-pos(3))/pathlength;
    
    % get sin & cos of azimuthal angle
    cp = 2;
    while ((cp >= 1) || (cp == 0));
        p = 2*rand - 1;
        q = 2*rand - 1;
        cp = q^2 + p^2;
    end
    cp = sqrt(cp);
    sinphi = q/cp;
    cosphi = p/cp;
    
    % get sin & cos of scatter angle
    costheta = 0.99*rand + 0.01;
    %costheta = 0.6173*rand + 0.38; % calculated for scatter angle of 67.5deg
    sintheta = sqrt(1-costheta^2);
    
    % change direction
    if abs(dir(3)) > 0.999 % if theta = 0, then equations simplify (also prevents error if costheta > 1)
        dir(1) = sintheta*cosphi;
        dir(2) = sintheta*sinphi;
        dir(3) = costheta*dir(3)/abs(dir(3));
    else
        sinnorm = sqrt(1-dir(3)^2);
        dirtemp1 = dir(1)*costheta + sintheta*(dir(1)*dir(3)*cosphi - dir(2)*sinphi)/sinnorm;
        dirtemp2 = dir(2)*costheta + sintheta*(dir(2)*dir(3)*cosphi - dir(1)*sinphi)/sinnorm;
        
        dir(3) = dir(3)*costheta - sinnorm*sintheta*cosphi;
        dir(1) = dirtemp1;
        dir(2) = dirtemp2;
        
        dnorm = sqrt(dir(1)^2+dir(2)^2+dir(3)^2);
        dir = dir./dnorm;
    end

end

%% Map how photon moves while in tissue
function [pos, dir, ftnwt, tissuesphere, nrus] = tissue(pos, dir, ftnwt, tissuesphere, nrus)
global nrel_at nrel_ta reflco portrad mut albedo stotal

    pathlength = -log(rand)/mut; %sample pathlength
    pos = adjust_pos(pos, dir, pathlength); %move to new position (adjust_pos)
    stotal = stotal + pathlength; % add scatter pathlength to total pathlength
    
    if pos(3) > 0 % then photon crossed surface
        pathlength = -pos(3)/dir(3); % backtrack to surface
        pos = adjust_pos(pos, dir, pathlength);
        
        if sqrt(pos(1)^2 + pos(2)^2) < portrad % then crossed inside port
            % reflect or refract
            bscwt = fresnel(nrel_ta, dir);
            
            if rand <= bscwt % reflects off back into tissue
                dir(3) = -dir(3);
            else % successfully transmits into sphere
%                 disp('Photon reentered sphere')
                dir = refract(nrel_at, dir);
                tissuesphere = 1;
            end
            
        else % crossed outside port
            % reflect or refract
            nrel = nrel_ta;
            bscwt = fresnel(nrel, dir);
            
            if rand <= bscwt % reflects off back into tissue
                dir(3) = -dir(3);
            else
                if sqrt(pos(1)^2 + pos(2)^2) < 15 % refraction occurs within cube surface area
                    ftnwt = ftnwt*reflco;
                    dir(3) = -dir(3);
                else % refraction occurs outside cube, kill photon
                    ftnwt = 0;
                end
            end
        end
    else % still in tissue
        ftnwt = ftnwt*albedo; % adjust weight
        dir = tis_scatter(dir); % sample new direction vector (tis_scatter)
        nrus = nrus + 1; % increment nrus
    end

end

%% Sample for new scatter angle in tissue
function dir = tis_scatter(dir)
global g g2

    % get sin & cos of azimuthal angle
    cp = 2;
    while ((cp >= 1) || (cp == 0));
        p = 2*rand - 1;
        q = 2*rand - 1;
        cp = q^2 + p^2;
    end
    cp = sqrt(cp);
    sinphi = q/cp;
    cosphi = p/cp;
    
    % get scatter angle from Henyey-Greenstein
    p = 2*g*rand + 1 - g;
    p = g2/p;
    p = p^2;
    
    costheta = (2 - g2 - p)/(2*g);
    sintheta = sqrt(1 - costheta^2);
    
    % change direction
    if abs(dir(3)) > 0.999 % if theta = 0, then equations simplify (also prevents error if costheta > 1)
        dir(1) = sintheta*cosphi;
        dir(2) = sintheta*sinphi;
        dir(3) = costheta*dir(3)/abs(dir(3));
    else
        sinnorm = sqrt(1-dir(3)^2);
        dirtemp1 = dir(1)*costheta + sintheta*(dir(1)*dir(3)*cosphi - dir(2)*sinphi)/sinnorm;
        dirtemp2 = dir(2)*costheta + sintheta*(dir(2)*dir(3)*cosphi - dir(1)*sinphi)/sinnorm;
        
        dir(3) = dir(3)*costheta - sinnorm*sintheta*cosphi;
        dir(1) = dirtemp1;
        dir(2) = dirtemp2;
        
        dnorm = sqrt(dir(1)^2+dir(2)^2+dir(3)^2);
        dir = dir./dnorm;
        
    end
    
end
