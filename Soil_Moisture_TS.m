%--------------------------BEGIN NOTE------------------------------%
% University of Virginia
%--------------------------END NOTE--------------------------------%
% INPUT DATA:
% All soil moisture data with 3d format (n:lat, m:lon, z: value)
%
% DESCRIPTION:
%
% REVISION HISTORY: 
% 9 Jul 2020 Hyunglok Kim; initial specification
%-----------------------------------------------------------------%
clc; clear

addpath('/sfs/qumulo/qproject/hydrosense/shape_files/sudan/') % shpae file for sudan
addpath('/sfs/qumulo/qproject/hydrosense/shape_files/moza/') % shape file for moza
shp_file_dir='/sfs/qumulo/qproject/hydrosense/shape_files/';

ifp='/project/hydrosense/matlab/mat/resampled_01/flood_movie/'; % input folder
ofp='/project/hydrosense/to_OUTSIDE/'; % output folder

%%
study_area='sudan';
no_data_mask=[0 0 0 0 0]; % if no data, set element with 1. [ASCAT, AMSR2, SMOS, SMAP, CYGNSS]
% As Venkat recommended, if you don't need CYGNSS, set no_data_mask as [0 0 0 0 1]
mv_day=2; %run_hour=3;

if strfind(study_area,'moza')
    year_=2019;
    start_doy=60; end_doy=95; % expected dates of flooding +/-
    mapping_start_doy=63; mapping_end_doy=85;
    a1=26;a2=30;
    b1=27;b2=34;
    r1=1; r2=1; title_position=0.85;
    
    iifp=[ifp,study_area,'_',num2str(year_),'/'];
    Run_domain_lower_left_lat=-22;
    Run_domain_lower_left_lon=32;
    Run_domain_upper_right_lat=-17;
    Run_domain_upper_right_lon=37;
    shp_file_name='Bairros.shp';
    
elseif strfind(study_area, 'sudan')
    year_=2019;
    start_doy=210; end_doy=250; % expected dates of flooding +/-
    mapping_start_doy=220; mapping_end_doy=244;
    
    % h-box scale
    a1=20;a2=55;
    b1=15; b2=45;
    r1=1.2; r2=2; title_position=0.85; %0.82 85
    
    %Khartoum City Flooding: (focused on the urban flooding)
    box1=[[32.4632,15.6186];[32.4632,15.4634]; [32.6243,15.4634]; [32.6242,15.6186]];
    % box scale    
    
    %%Khartoum Greater Area Flooding: (focused on all flooding including agriculture in south)
    box2=[[32.0517,15.6254];[32.0517,13.6858];[33.5596,13.6858];[33.5596,15.6254]];    

    iifp=[ifp,study_area,'_',num2str(year_),'/'];
    Run_domain_lower_left_lat=12;
    Run_domain_lower_left_lon=30;
    Run_domain_upper_right_lat=18;
    Run_domain_upper_right_lon=35;
    shp_file_name='sdn_admbndl_admALL_cbs_nic_ssa_itos_20200317.shp';
end

hh_start_doy=(start_doy-1)*48+1; %half-hour temporal resolution
hh_end_doy=end_doy*48;

% make title
title_list=[]; k=1;
title_num=[];
for i=hh_start_doy:hh_end_doy
    
    t=datevec(datenum(year_,1,1)+i/48);
    t_year=num2str(t(1));
    t_month=num2str(t(2));
    if t(2)<10
        t_month=['0',num2str(t(2))];
    end
    
    t_day=num2str(t(3));
    if t(3)<10
        t_day=['0',num2str(t(3))];
    end
    
    t_hr=num2str(t(4));
    if t(4)<10
        t_hr=['0',num2str(t(4))];
    end
    
    t_mn=num2str(t(5));
    if t(5)<10
        t_mn=['0',num2str(t(5))];
    end
    title_list{k,1}=[t_year,'/',t_month,'/',t_day,'/',t_hr,':',t_mn];
    title_num(k,1)=datenum(t);
    k=k+1;
    
end

% load all data
t_SM=load([iifp,'ASCAT/ASCAT_SM_',num2str(year_),'.mat']);
ASCAT_SM=t_SM.ASCAT_SM;

t_SM=load([iifp,'AMSR2/AMSR2_SM_',num2str(year_),'_X_DES.mat']);
AMSR2_SM_DES=t_SM.AMSR2_SM;

t_SM=load([iifp,'AMSR2/AMSR2_SM_',num2str(year_),'_X_ASC.mat']);
AMSR2_SM_ASC=t_SM.AMSR2_SM;

t_SM=load([iifp,'SMOS/SMOS_SM_',num2str(year_),'_DES.mat']);
SMOS_SM_DES=t_SM.SMOS_SM;

t_SM=load([iifp,'SMOS/SMOS_SM_',num2str(year_),'_ASC.mat']);
SMOS_SM_ASC=t_SM.SMOS_SM;

t_SM=load([iifp,'SMAP/SMAP_SM_',num2str(year_),'_DES.mat']);
SMAP_SM_DES=t_SM.SMAP_SM;

t_SM=load([iifp,'SMAP/SMAP_SM_',num2str(year_),'_ASC.mat']);
SMAP_SM_ASC=t_SM.SMAP_SM;

t_SM=load([iifp,'CYGNSS/CYGNSS_SM_',num2str(year_),'_sub_daily.mat']);
CYGNSS_SM=t_SM.CYGNSS_SM;

% rescale to the reference data
ref_data=ASCAT_SM; % it can be any data

r_ASCAT_SM=[];
r_AMSR2_SM_DES=[]; r_AMSR2_SM_ASC=[];
r_SMOS_SM_DES=[]; r_SMOS_SM_ASC=[];
r_SMAP_SM_DES=[]; r_SMAP_SM_ASC=[];
r_CYGNSS_SM=[];

for i=1:size(ref_data,1)
    for j=1:size(ref_data,2)
        t_ref=squeeze(ref_data(i,j,:));
        
        t_sat=squeeze(ASCAT_SM(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_ASCAT_SM(i,j,:)=t_sat;
        else
            r_ASCAT_SM(i,j,:)=t_nor;
        end
        
        t_sat=squeeze(AMSR2_SM_DES(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_AMSR2_SM_DES(i,j,:)=t_sat;
        else
            r_AMSR2_SM_DES(i,j,:)=t_nor;
        end
        
        t_sat=squeeze(AMSR2_SM_ASC(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_AMSR2_SM_ASC(i,j,:)=t_sat;
        else
            r_AMSR2_SM_ASC(i,j,:)=t_nor;
        end
        
        t_sat=squeeze(SMOS_SM_DES(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_SMOS_SM_DES(i,j,:)=t_sat;
        else
            r_SMOS_SM_DES(i,j,:)=t_nor;
        end
        
        t_sat=squeeze(SMOS_SM_ASC(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_SMOS_SM_ASC(i,j,:)=t_sat;
        else
            r_SMOS_SM_ASC(i,j,:)=t_nor;
        end
        
        t_sat=squeeze(SMAP_SM_ASC(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_SMAP_SM_ASC(i,j,:)=t_sat;
        else
            r_SMAP_SM_ASC(i,j,:)=t_nor;
        end
                        
        t_sat=squeeze(SMAP_SM_DES(i,j,:));
        t_nor=SAT_NOR([t_ref, t_sat]);
        if sum(isnan(t_nor))==0
            r_SMAP_SM_DES(i,j,:)=t_sat;
        else
            r_SMAP_SM_DES(i,j,:)=t_nor;
        end
        
        t_sat=squeeze(CYGNSS_SM(i,j,:));
        t_nor=SAT_NOR_V2(t_ref, t_sat);
        if sum(isnan(t_nor))==0
            r_CYGNSS_SM(i,j,:)=t_sat;
        else
            r_CYGNSS_SM(i,j,:)=t_nor;
        end
    end
end

% Fill soil moisture data and make half-houlry data (basic)
hh_SM=nan([size(ref_data,1), size(ref_data,2), size(ref_data,3)*48]);

if no_data_mask(1)~=1
    hh_SM(:,:,44:48:end)=r_ASCAT_SM;
end

if no_data_mask(2)~=1
    hh_SM(:,:,4:48:end)=r_AMSR2_SM_DES;
    hh_SM(:,:,28:48:end)=r_AMSR2_SM_ASC;
end

rr_SMAP_SM_ASC=r_SMAP_SM_ASC;
rr_SMAP_SM_DES=r_SMAP_SM_DES;
if no_data_mask(3)~=1
    nan_data=isnan(rr_SMAP_SM_DES);
    rr_SMAP_SM_DES(nan_data)=r_SMOS_SM_ASC(nan_data);
    nan_data=isnan(rr_SMAP_SM_ASC);
    rr_SMAP_SM_ASC(nan_data)=r_SMOS_SM_DES(nan_data);
    
    hh_SM(:,:,13:48:end)=r_SMOS_SM_ASC;
    hh_SM(:,:,37:48:end)=r_SMOS_SM_DES;
end

if no_data_mask(4)~=1
    hh_SM(:,:,13:48:end)=rr_SMAP_SM_DES;
    hh_SM(:,:,37:48:end)=rr_SMAP_SM_ASC;
end

if no_data_mask(5)~=1
    hh_SM(:,:,7:48:end)=r_CYGNSS_SM(:,:,1:4:end);
    hh_SM(:,:,19:48:end)=r_CYGNSS_SM(:,:,2:4:end);
    hh_SM(:,:,31:48:end)=r_CYGNSS_SM(:,:,3:4:end);
    hh_SM(:,:,43:48:end)=r_CYGNSS_SM(:,:,4:4:end);
end

% fill with previous day data
f_SM=hh_SM;
for i=1:size(f_SM,3)-1
    t_SM=hh_SM(:,:,i+1);
    tt_SM=f_SM(:,:,i);
    nodata=isnan(t_SM);
    t_SM(nodata)=tt_SM(nodata);    
    f_SM(:,:,i+1)=t_SM;    
end

% calculate moving average data
[~, m_SM]=ano_coarse(f_SM(:,:,hh_start_doy:hh_end_doy), mv_day*48);

%% time series avg SM
% calculate 95% of SM
t_ASCAT_SM=ASCAT_SM(a1:a2, b1:b2,:);
foo_SM=[];
for i=1:size(t_ASCAT_SM,3)
    t=t_ASCAT_SM(:,:,i);
    foo_SM(i)=nanmean(t(:));
end
%t_ASCAT_SM=t_ASCAT_SM(:);
foo_SM=sort(foo_SM);
p_ASCAT_SM=[];
for i=1:numel(foo_SM)
    p_ASCAT_SM(i)=i/numel(foo_SM);
end
t_max=foo_SM(min((p_ASCAT_SM>0.95)));

t_avg=[]; k=1;
 for i=(mapping_start_doy-start_doy)*48+1:(mapping_end_doy-start_doy+1)*48
     t=(m_SM(a1:a2,b1:b2,i));
     t_avg(k)=nanmean(t(:));
     k=k+1;
 end
 %% plot the soil moisture TS
 plot(title_num((mapping_start_doy-start_doy)*48+1:(mapping_end_doy-start_doy+1)*48),t_avg, 'linewidth',1, 'color','k')
 ylim([0.2 0.35])
 datetick('x','mm-dd')
 set(gca,'fontweight','bold', 'fontsize',20)
 ylabel('Soil Moisture (m^3/m^3)','fontsize',30)
 ay=get(gca,'YtickLabel');
 set(gca, 'YtickLabel',ay,'fontsize',25)
 grid on
 
 set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
 title({'Khartoum (Urban+Agriculture)','-Soil Moisture'},'fontsize',35)
 %title('Khartoum (Urban)-Soil Moisture','fontsize',35)
 