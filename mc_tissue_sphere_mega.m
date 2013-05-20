%% Monte Carlo simulation
%for measuring photons detected with sphere over tissue
%used to calculate scattering losses as a function of optical properties

%% MAIN PROGRAM
function [] = mc_tissue_sphere_mega

params = struct; %Defines empty struct
% seed the random number generator based on the current time
%rng('shuffle'); %Works for new matlabs

stream = RandStream('mt19937ar','Seed',sum(100*clock));  %Needed for older matlabs
RandStream.setDefaultStream(stream); %Needed for older matlabs

% read XLSX file containing list of MC sims to run & related data
[file, nruns] = read_list_sims;
disp(['Processing file ',file])
% for each simulation,

for currun = 1:nruns %first run starts at 1 (Sheet1)

    params = read_param(currun, file); %Fill up parameter struct
    stotalbin = zeros(10001,1); %Zero statistics for this run
    disp(['Currently processing run ', int2str(currun), '/', int2str(nruns), ' ', int2str(1000*params.kftn), ' Photons'])
    tic %Start timer for performance measurements
    parfor_progress(1000*params.kftn); %Terminal-based progress indicator
        parfor ftncount = 1:1000*params.kftn
        
            local_stotalbin = zeros(10001,1); %Initialize statistics bin for *just* this photon
            [pos, dir, nrus, ftnwt, tissuesphere, stotal] = ftnini(params); % initialize photon, these are the per-photon parameters
            
            while (ftnwt > 1e-16)
                if tissuesphere == 1
                    [pos, dir, nrus, ftnwt, tissuesphere, stotal, local_stotalbin] = sphere(pos, dir, nrus, ftnwt, tissuesphere, stotal, local_stotalbin, params);
                elseif tissuesphere == -1
                    [pos, dir, nrus, ftnwt, tissuesphere, stotal] = tissue(pos, dir, nrus, ftnwt, tissuesphere, stotal, params);
                else
                    disp('Something has gone terribly wrong (tissuesphere =/ +/-1).')
                    ftnwt = 0;
                end
                
                if (pos(3) < 0) && (norm(pos) >= 100) % too far away
                    ftnwt = 0;
                    %radius = norm(pos)
                    %disp('Too far away')
                end
                
                %if nrus >= params.rusnum % too many scatters
                %    [ftnwt,nrus] = rusrul(ftnwt,nrus,params);
                %end
            end
            stotalbin = stotalbin + local_stotalbin; %Accumulate results into main statistics for run
            parfor_progress; %Update progress bar
        end
        parfor_progress(0); %Cleanup progress bar files
        toc %Stop timer for performance measurements, output time
        str = num2str(currun);
        runloc = strcat('Sheet',str);
        xlswrite(file,stotalbin,runloc,'D1')% print detectwt to current run's sheet
 
    
end

end

%% Read in location of data and number of simulations to run
function [file, nruns] = read_list_sims % VERIFIED

[filename, pathname] = uigetfile('*.xlsx', 'Select the XLSX file containing the MC data to run');
file = strcat(pathname,filename);
nruns = xlsread(file,'Sheet0');

end

%% Read in parameters for the current simulation
function params = read_param(currun, file) % VERIFIED

    str = num2str(currun);
    runloc = strcat('Sheet',str);
    [num_data, text_data] = xlsread(file,runloc);
    
    params.g = num_data(1); % anisotropy parameter
    params.g2 = 1- params.g^2;
    params.nrel_at = num_data(2); % relative index of refraction (air to tissue)
    params.nrel_ta = 1/params.nrel_at;
    params.reflco = num_data(3); % reflection coefficient of sphere
    params.rusnum = num_data(4); % # of scatters before Russian roulette
    params.rusfrac = num_data(5); % survival fraction applied to successful roulettes
    params.kftn = num_data(6); % number of groups of thousands of photons
    params.israd = num_data(7); % radius of integrating sphere (mm)
    params.portrad = num_data(8); % radius of detection port (mm)
    params.midz = sqrt(params.israd^2-params.portrad^2); % height of center of sphere
    params.coshalf = num_data(9); % cosine of the fiber half angle
    params.surfrad = (sqrt(1-params.coshalf^2)/params.coshalf)*(params.israd + params.midz); % tissue surface "seen" by detector fiber
    params.mua = num_data(10); % absorption coefficient (mm^-1)
    params.musp = num_data(11); % reduced scatter coefficient (mm^-1)
    params.mut = params.mua + params.musp/(1-params.g);
    params.albedo = (params.mut-params.mua)/params.mut; % transport/scatter albedo
    
end

%% Initialize photon
function [pos, dir, nrus, ftnwt, tissuesphere, stotal] = ftnini(params) % VERIFIED

pos = [params.israd 0 params.midz];

costheta = params.coshalf + (1 - params.coshalf)*rand;
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
function [ftnwt,nrus] = rusrul(ftnwt, nrus, params)

    if (rand <= params.rusfract)
        nrus = -10*nrus;
        ftnwt = ftnwt/params.rusfrac;
    else
        ftnwt = 0;
        disp('Roulette killed.')
    end
    
end

%% Map how photon moves while in sphere
function [pos, dir, nrus, ftnwt, tissuesphere, stotal, stotalbin] = sphere(pos, dir, nrus, ftnwt, tissuesphere, stotal, stotalbin, params)

    oldpos = pos;
    pos = raytrace(pos, dir, params); % raytrace to new position on sphere
    
    if pos(3) <= 0 % then have gone into tissue
        % backtrack to tissue surface
        pathlength = -pos(3)/dir(3);
        pos = adjust_pos(pos, dir, pathlength);
        
        % at tissue surface, check to see if photon reflects or transmits:
        bscwt = fresnel(params.nrel_at, dir);
        
        if rand <= bscwt % reflects off surface
            dir(3) = -dir(3);
        else % successfully transmits into tissue
            dir = refract(params.nrel_ta, dir);
            tissuesphere = -1;
        end
        
    else % photon is still in sphere
        if (sqrt(oldpos(1)^2+oldpos(2)^2) < params.surfrad) && (abs(oldpos(3)) < 0.1) && (pos(3) > (params.midz + sqrt(params.israd^2 - 2.5^2)))
%  disp('Photon scored')            
            if stotal <= 0
                stotalbin(1) = stotalbin(1) + ftnwt;
            elseif stotal > 1000
				ftnwt = 0;
                %stotalbin(10001) = stotalbin(10001) + ftnwt;
            else
                bin = ceil(stotal*10 + 1); % find bin
                stotalbin(bin) = stotalbin(bin) + ftnwt; % add total pathlength of photon to bin
            end
            ftnwt = 0;
            
        else % is reflected off sphere wall
            ftnwt = ftnwt*params.reflco; % adjust photon weight
            dir = sphere_scat(pos, params); % sample for new direction vector
        end
    end

end

%% Calculates the pathlength and then position along the current direction vector to the next spot on the sphere wall
function newpos = raytrace(pos, dir, params) % VERIFIED

    a = 1;
    b = 2*(pos(1)*dir(1) + pos(2)*dir(2) + (pos(3) - params.midz)*dir(3));
    c = pos(1)^2 + pos(2)^2 + (pos(3) - params.midz)^2 - params.israd^2;
    
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
function dir = sphere_scat(pos, params)

    % point toward center of sphere
    pathlength = sqrt(pos(1)^2 + pos(2)^2 + (params.midz-pos(3))^2);
    dir(1) = -pos(1)/pathlength;
    dir(2) = -pos(2)/pathlength;
    dir(3) = (params.midz-pos(3))/pathlength;
    
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
        
        dir = dir./norm(dir);
    end

end

%% Map how photon moves while in tissue
function [pos, dir, nrus, ftnwt, tissuesphere, stotal] = tissue(pos, dir, nrus, ftnwt, tissuesphere, stotal, params)

    pathlength = -log(rand)/params.mut; %sample pathlength
    pos = adjust_pos(pos, dir, pathlength); %move to new position (adjust_pos)
    stotal = stotal + pathlength; % add scatter pathlength to total pathlength
    
    if pos(3) > 0 % then photon crossed surface
        pathlength = -pos(3)/dir(3); % backtrack to surface
		stotal = stotal - pathlength % correct pathlength by subtracting distance travelled in air
        pos = adjust_pos(pos, dir, pathlength);
        
        if sqrt(pos(1)^2 + pos(2)^2) < params.portrad % then crossed inside port
            % reflect or refract
            bscwt = fresnel(params.nrel_ta, dir);
            
            if rand <= bscwt % reflects off back into tissue
                dir(3) = -dir(3);
            else % successfully transmits into sphere
%                 disp('Photon reentered sphere')
                dir = refract(params.nrel_at, dir);
                tissuesphere = 1;
            end
            
        else % crossed outside port
            % reflect or refract

            bscwt = fresnel(params.nrel_ta, dir);
            
            if rand <= bscwt % reflects off back into tissue
                dir(3) = -dir(3);
            else
                if sqrt(pos(1)^2 + pos(2)^2) < 25 % refraction occurs within cube surface area
                    ftnwt = ftnwt*params.reflco;
					costheta = -rand;
					sintheta = sqrt(1-costheta^2);
					
					cp = 2;
					while ((cp >= 1) || (cp == 0));
						p = 2*rand - 1;
						q = 2*rand - 1;
						cp = q^2 + p^2;
					end
					cp = sqrt(cp);
					sinphi = q/cp;
					cosphi = p/cp;
					
					dir(1) = sintheta*cosphi;
					dir(2) = sintheta*sinphi;
					dir(3) = costheta;
                    %dir(3) = -dir(3);
                else % refraction occurs outside cube, reflect but modify weight by backscatter factor
                    dir(3) = -dir(3);
					ftnwt = ftnwt*bscwt;
                end
            end
        end
    else % still in tissue
        ftnwt = ftnwt*params.albedo; % adjust weight
        dir = tis_scatter(dir, params); % sample new direction vector (tis_scatter)
        nrus = nrus + 1; % increment nrus
    end

end

%% Sample for new scatter angle in tissue
function dir = tis_scatter(dir, params)

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
    p = 2*params.g*rand + 1 - params.g;
    p = params.g2/p;
    p = p^2;
    
    costheta = (2 - params.g2 - p)/(2*params.g);
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
        
        dir = dir./norm(dir);
    end
    
end
