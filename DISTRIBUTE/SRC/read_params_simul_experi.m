function [experiment,experiment_hadc,experiment_btpde] ...
    = read_params_simul_experi(fname_experiment)

% read experiment parameters
% 
% Input:
%         fname_experiment
%             
% Output:
% 
%     1. experiment is a structure with 6 elements:
%         ngdir_total 
%         gdir        
%         sdeltavec   
%         bdeltavec    
%         seqvec       
%         npervec      
%         
%     2. experiment_hadc is a structure with 8 elements:
%         ngdir_total 
%         gdir         
%         sdeltavec    
%         bdeltavec    
%         seqvec       
%         npervec     
%         rtol        
%         atol        
% 
%     3. experiment_btpde is a structure with 10 elements:
%         ngdir_total 
%         gdir         
%         sdeltavec    
%         bdeltavec    
%         seqvec       
%         npervec     
%         rtol        
%         atol        
%         qvalues     
%         bvalues     

SEQ_DEFINITIONS
global QVAL
global BVAL
global BDELTA SDELTA SEQ OGSEPER 

ndim = 3;
clear bvec;

fid=fopen(fname_experiment);

tline = fgetl(fid);
experiment.ngdir_total= sscanf(tline,'%f',1);

tline = fgetl(fid);
experiment.gdir= sscanf(tline,'%f',ndim);
experiment.gdir = experiment.gdir/norm(experiment.gdir);

tline = fgetl(fid);
nd= sscanf(tline,'%f',1);

tline = fgetl(fid);
experiment.sdeltavec= sscanf(tline,'%f',nd);

tline = fgetl(fid);
experiment.bdeltavec= sscanf(tline,'%f',nd);

tline = fgetl(fid);
experiment.seqvec= sscanf(tline,'%f',nd);
tline = fgetl(fid);
experiment.npervec= sscanf(tline,'%f',nd);

tline = fgetl(fid);
do_hadc = sscanf(tline,'%f',1);
tline = fgetl(fid);
atmp= sscanf(tline,'%f',2);
if (do_hadc ~= 0)
	experiment_hadc = experiment;
	experiment_hadc.rtol= atmp(1);
	experiment_hadc.atol = atmp(2);
else
	experiment_hadc = [];
end

tline = fgetl(fid);
do_btpde = sscanf(tline,'%f',1);
tline = fgetl(fid);
atmp= sscanf(tline,'%f',2);

if (do_btpde ~= 0)
	experiment_btpde = experiment;
	experiment_btpde.rtol = atmp(1);
	experiment_btpde.atol = atmp(2);
	tline = fgetl(fid);
	nb = sscanf(tline,'%f',1);

	tline = fgetl(fid);
	use_blimits= sscanf(tline,'%f',1);

	tline = fgetl(fid);
	const_q = sscanf(tline,'%f',1);

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
else
	experiment_btpde = [];
end

fclose(fid);

if (do_btpde ~= 0)
	nb = max(length(qvec),length(bvec));
	nexperi = length(experiment_btpde.sdeltavec);

	experiment_btpde.qvalues = zeros(nexperi,nb);
	experiment_btpde.bvalues = zeros(nexperi,nb);

	for iexperi = 1:nexperi
		SDELTA = experiment_btpde.sdeltavec(iexperi);
		BDELTA = experiment_btpde.bdeltavec(iexperi);
		TE = SDELTA+BDELTA;
		SEQ = experiment_btpde.seqvec(iexperi); % for choosing case PGSE, OGSEcos or OGSEsin
		omega = 2*pi*experiment_btpde.npervec(iexperi)/SDELTA;
		OGSEPER = 1./omega*2*pi;%% set up number for OGSE
		
		for ib = 1:nb
			if (length(bvec) == 0)
				experiment_btpde.qvalues(iexperi,ib) = qvec(ib);
				QVAL = experiment_btpde.qvalues(iexperi,ib);
				
				experiment_btpde.bvalues(iexperi,ib) = seqbvaluenoq*experiment_btpde.qvalues(iexperi,ib)^2;
			else
				if (iexperi == 1 | const_q ~= 1)
					experiment_btpde.bvalues(iexperi,ib) = bvec(ib);
					BVAL = experiment_btpde.bvalues(iexperi,ib);
					experiment_btpde.qvalues(iexperi,ib)  = sqrt(experiment_btpde.bvalues(iexperi,ib)/seqbvaluenoq); % water proton gyromagnetic ratio* norm of g
				else
					experiment_btpde.qvalues(iexperi,ib) = experiment_btpde.qvalues(1,ib);
					QVAL = experiment_btpde.qvalues(iexperi,ib);
					experiment_btpde.bvalues(iexperi,ib) = seqbvaluenoq*experiment_btpde.qvalues(iexperi,ib)^2;
				end
			end
		end
    end
end

 