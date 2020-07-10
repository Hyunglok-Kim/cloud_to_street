function SATNOR=SAT_NOR_V2(ref, target)
%--------------------------BEGIN NOTE------------------------------%
% University of Virginia
%--------------------------END NOTE--------------------------------%
% ARGUMENTS:
%
% DESCRIPTION:
% Normalize target data into reference data
% nan value can be included
%
% REVISION HISTORY: 
% 9 Jul 2020 Hyunglok Kim; initial specification
%-----------------------------------------------------------------%

nod_1=sum(~isnan(ref)); % number of non-nan value
nod_2=sum(~isnan(target)); % number of non-nan value

if nod_1 < 15 || nod_2 < 15
    disp([num2str(nod_1),'/',num2str(nod_2) '<-number of data is too small']);
    SATNOR=nan(size(target));
    return
end

OBS=ref;SAT=target;
OBS_std=nanstd(OBS);
OBS_mean=nanmean(OBS);
SAT_std=nanstd(SAT);
SAT_mean=nanmean(SAT);

SATNOR=(SAT-SAT_mean)*OBS_std/SAT_std+OBS_mean;

    