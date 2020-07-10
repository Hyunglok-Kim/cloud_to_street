function [ano, coarse]=ano_coarse(dataset, moving_day)
%--------------------------BEGIN NOTE------------------------------%
% University of Virginia
%--------------------------END NOTE--------------------------------%
% ARGUMENTS:
% dataset: n x m x val 
% moving_day: target moving average day 
% DESCRIPTION:
% Normalize target data into reference data
% Based on Su, C. H.,et al. (2014).
%
% REVISION HISTORY: 
% 6 Apr 2016 Hyunglok Kim; initial specification
%-----------------------------------------------------------------%

nod=nan(size(dataset));
coarse=nod;ano=nod;
for doy=(floor(moving_day/2)+1):size(dataset,3)-floor(moving_day/2)

    a1=doy-floor(moving_day/2);
    a2=doy+floor(moving_day/2);
    temp_1=dataset(:,:,a1:a2);
    temp_2=sum(~isnan(temp_1),3);

    temp_3=repmat(temp_2<5, [1 1 size(temp_1,3)]);
    temp_1(temp_3)=nan;

    temp_2(temp_2<5)=nan;
    nod(:,:,doy)=temp_2;
    temp_4=nansum(temp_1,3);

    coarse(:,:,doy)=temp_4./(temp_2+1);
end
ano=dataset-coarse;
%coarse(isnan(ano))=nan;
