function [DEFF] = deff_sta(D0,VOL,SAu,SDELTA,BDELTA,seq,nper)

global PGSE OGSEsin OGSEcos dPGSE
global SEQ OGSEPER 


if (seq == PGSE)
    DEFF = D0-D0^(3/2)*16/35/VOL/sqrt(pi)/SDELTA^2/(3*BDELTA-SDELTA)*...
        (-2*(SDELTA^(7/2)+BDELTA^(7/2))+(BDELTA-SDELTA)^(7/2)...
        +(BDELTA+SDELTA)^(7/2))*SAu;
else
	SEQ = seq	
    omega = 2*pi*nper/SDELTA;
    OGSEPER = 1./omega*2*pi;%% set up number for OGSE
	difftime = seqdifftime
    DEFF = D0*(1-4/3/sqrt(pi)*sqrt(D0)*SAu/VOL*sqrt(difftime));
end
