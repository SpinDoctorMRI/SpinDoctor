function [gdir,bvalues,qvalues,sdeltavec,bdeltavec,seqvec,npervec,...
    rtol_bt,atol_bt,rtol_deff,atol_deff,const_q,tetgen_cmd] ...
    = read_simulation_parameters_experiment(fname_experiment);

global QVAL
global BVAL
global BDELTA SDELTA SEQ OGSEPER 

ndim = 3;
clear bvec;

fid=fopen(fname_experiment);

tline = fgetl(fid);
gdir= sscanf(tline,'%f',ndim);

tline = fgetl(fid);
nb= sscanf(tline,'%f',1);

tline = fgetl(fid);
use_blimits= sscanf(tline,'%f',1);

if (use_blimits == 0)
    tline = fgetl(fid);
    bvec= sscanf(tline,'%f',nb);
    qvec = [];
elseif (use_blimits == 1)
    tline = fgetl(fid);
    blimits = sscanf(tline,'%f',2);
    bvec = linspace(blimits(1),blimits(2),nb);
    qvec = [];
elseif (use_blimits == 2)
    tline = fgetl(fid);
    glimits = sscanf(tline,'%f',2);
    qvec = linspace(glimits(1),glimits(2),nb);
    bvec = [];
end    

tline = fgetl(fid);
nd= sscanf(tline,'%f',1);

tline = fgetl(fid);
sdeltavec= sscanf(tline,'%f',nd);

tline = fgetl(fid);
bdeltavec= sscanf(tline,'%f',nd);

tline = fgetl(fid);
seqvec= sscanf(tline,'%f',nd);
tline = fgetl(fid);
npervec= sscanf(tline,'%f',nd);

tline = fgetl(fid);
atmp= sscanf(tline,'%f',2);
rtol_bt= atmp(1);
atol_bt = atmp(2);

tline = fgetl(fid);
atmp= sscanf(tline,'%f',2);
rtol_deff= atmp(1);
atol_deff = atmp(2);


tline = fgetl(fid);
const_q = sscanf(tline,'%f',1);

tline = fgetl(fid);
[strpos] = regexp(tline,"'");
tetgen_cmd = tline(strpos(1)+1:strpos(2)-1);

fclose(fid);

gdir = gdir/norm(gdir);

nb = max(length(qvec),length(bvec));
nexperi = length(sdeltavec);

qvalues = zeros(nexperi,nb);
bvalues = zeros(nexperi,nb);

for iexperi = 1:nexperi
    SDELTA = sdeltavec(iexperi);
    BDELTA = bdeltavec(iexperi);
    TE = SDELTA+BDELTA;
    SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
    omega = 2*pi*npervec(iexperi)/SDELTA;
    OGSEPER = 1./omega*2*pi;%% set up number for OGSE
    
    for ib = 1:nb
        if (length(bvec) == 0)
            qvalues(iexperi,ib) = qvec(ib);
            QVAL = qvalues(iexperi,ib);
            
            bvalues(iexperi,ib) = seqbvaluenoq*qvalues(iexperi,ib)^2;
        else
            if (iexperi == 1 | const_q ~= 1)
                bvalues(iexperi,ib) = bvec(ib);
                BVAL = bvalues(iexperi,ib);
                qvalues(iexperi,ib)  = sqrt(bvalues(iexperi,ib)/seqbvaluenoq); % water proton gyromagnetic ratio* norm of g
            else
                qvalues(iexperi,ib) = qvalues(1,ib);
                QVAL = qvalues(iexperi,ib);
                bvalues(iexperi,ib) = seqbvaluenoq*qvalues(iexperi,ib)^2;
            end
        end
    end
end


 