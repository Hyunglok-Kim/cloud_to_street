function SATNOR=SAT_NOR(OBS_SAT)
%--------------------------BEGIN NOTE------------------------------%
% University of Virginia
%--------------------------END NOTE--------------------------------%
% ARGUMENTS:
% OBS_SAT = N X 2 array 
%
% DESCRIPTION:
% Normalize target data into reference data
% nan value can be included
%
% REVISION HISTORY: 
% 9 Jul 2020 Hyunglok Kim; initial specification
%-----------------------------------------------------------------%

if size(OBS_SAT,2)~=2
    disp(['must composed with N X 2 array'])
    SATNOR='bad data array';
    return
end

nod=sum(~isnan(sum(OBS_SAT,2))); % number of non-nan value
if nod < 15
    disp([num2str(nod),'<-number of data is too small']);
    SATNOR=nan(size(OBS_SAT(:,1)));
    return
end

OBS=OBS_SAT(:,1);
SAT=OBS_SAT(:,2);
OBS_std=nanstd(OBS);
OBS_mean=nanmean(OBS);
SAT_std=nanstd(SAT);
SAT_mean=nanmean(SAT);

SATNOR=(SAT-SAT_mean)*OBS_std/SAT_std+OBS_mean;

    