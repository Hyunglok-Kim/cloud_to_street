clc; clear all

addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/SAT_data_related_CODE/');
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/mapping_code/');
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs')
addpath('/sfs/qumulo/qproject/hydrosense/shape_files/sudan/') % shpae file of sudan
addpath('/sfs/qumulo/qproject/hydrosense/shape_files/moza/') % shape file of moza
shp_file_dir='/sfs/qumulo/qproject/hydrosense/shape_files/';
ifp='/project/hydrosense/matlab/mat/resampled_01/flood_movie/'; % input folder
ofp='/project/hydrosense/to_OUTSIDE/'; % output folder

%%
study_area='sudan';
no_data_mask=[0 0 0 0 0]; % if no data, set to 1 ASCAT, AMSR2, SMOS, SMAP, CYGNSS
mv_day=2;
run_hour=3;

if strfind(study_area,'moza')
    year_=2019;
    start_doy=60; end_doy=95;
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
    start_doy=210; end_doy=250; %210; 250
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
Run_domain_resolution_dx=0.1;
Run_domain_resolution_dy=0.1;

hh_start_doy=(start_doy-1)*48+1;
hh_end_doy=end_doy*48;

temp_lat=Run_domain_upper_right_lat:-Run_domain_resolution_dy:Run_domain_lower_left_lat;
temp_lon=Run_domain_lower_left_lon:Run_domain_resolution_dx:Run_domain_upper_right_lon;
[target_lon,target_lat]=meshgrid(temp_lon,temp_lat);

% make title
title_list=[]; k=1;
title_num=[];
for i=hh_start_doy:hh_end_doy
    %for j=1:48
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
    %end
end
% load data
clc
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

% rescale to reference data. 
ref_data=ASCAT_SM;%it can be any data; SMAP_SM_DES;

r_ASCAT_SM=[];
r_AMSR2_SM_DES=[]; r_AMSR2_SM_ASC=[];
r_SMOS_SM_DES=[]; r_SMOS_SM_ASC=[];
r_SMAP_SM_DES=[]; r_SMAP_SM_ASC=[];
r_CYGNSS_SM=[];
%
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

% mkae half-houlry data
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
%
% t=squeeze(hh_SM(10,10,:));
% plot(1:numel(t), t,'o');
% hold on
%t=squeeze(m_SM(10,10,:));
%plot(1:numel(t), t,'x');
% check spatail averaged SM data
%% make movie
%a1=24;a2=27;
%b1=25;b2=28;

a1=23;a2=45;
b1=20;b2=38;

i=(mapping_start_doy-start_doy)*48+1; %mv_day/2*48+1;
target=m_SM(a1:a2,b1:b2,i);
%max_t=max(target(:));
%min_t=min(target(:));

target(target>=0.4)=0.4; target(target<0)=0;
%target(1)=min_t; target(2)=max_t;
target(1)=0; target(2)=0.4;
%%
Statistic_Mapping_10km_flood(target,target_lat(a1:a2,b1:b2),target_lon(a1:a2,b1:b2));
load coastlines
plotm(coastlat, coastlon, 'linewidth',1,'color','k')
title(title_list{i,1},'Units','normalized','position',[0.5 title_position],'fontsize',20, 'fontname','times new roman')

t=shaperead([shp_file_dir,study_area,'/',shp_file_name]);
map_x=[]; map_y=[];

for i=1:size(t)
    t_x=t(i).X; t_y=t(i).Y;
    map_x=[map_x, t_x];
    map_y=[map_y, t_y];
    hold on
end
plotm(map_y, map_x,'color','k')
hold on

for i=1:4
    geoshow(box2(i,2),box2(i,1),'DisplayType','Point','Marker','p','markeredgecolor',[0 0 0],'markersize',15,'linewidth',2.5)
    geoshow(box1(i,2),box1(i,1),'DisplayType','Point','Marker','.','markeredgecolor','black','markersize',15,'linewidth',2.5)
end
%
set(gca, 'plotboxaspectratio',[r1,r2,1])
 
%geoshow(15.25,32.69,'DisplayType','Point','Marker','p','markeredgecolor','black','markersize',25,...
%    'linewidth',2.5)

gif([ofp, study_area,'_SM_city_',num2str(no_data_mask(1)),num2str(no_data_mask(2)),num2str(no_data_mask(3)),...
    num2str(no_data_mask(4)),num2str(no_data_mask(5)),'.gif'],'DelayTIme',0.1)
for i=(mapping_start_doy-start_doy)*48+2:run_hour*2:(mapping_end_doy-start_doy+1)*48
    target=m_SM(a1:a2,b1:b2,i);
    target(target>0.4)=0.4;
    surfm(target_lat(a1:a2,b1:b2), target_lon(a1:a2,b1:b2), target)
    %geoshow(15.25,32.69,'DisplayType','Point','Marker','p','markeredgecolor','black','markersize',25,...
    %    'linewidth',2.5)
    title(title_list{i,1},'Units','normalized','position',[0.5 title_position 1],'fontsize',20, 'fontname','times new roman')
    %load coastlines
    %plotm(coastlat, coastlon, 'linewidth',1.5,'color','k')
    plotm(map_y, map_x,'color','k')
    %set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %set(h, 'title',title_list{i,1})
    
    for i=1:4
        geoshow(box2(i,2),box2(i,1),'DisplayType','Point','Marker','p','markeredgecolor',[0 0 0],'markersize',15,'linewidth',2.5)
        geoshow(box1(i,2),box1(i,1),'DisplayType','Point','Marker','.','markeredgecolor','black','markersize',15,'linewidth',2.5)
    end

    i
    gif
end
%% time series avg SM
% calculate 95% of SM
t_ASCAT_SM=ASCAT_SM(a1:a2, b1:b2,:);
tt__SM=[];
for i=1:size(t_ASCAT_SM,3)
    t=t_ASCAT_SM(:,:,i);
    tt__SM(i)=nanmean(t(:));
end
%t_ASCAT_SM=t_ASCAT_SM(:);
tt__SM=sort(tt__SM);
p_ASCAT_SM=[];
for i=1:numel(tt__SM)
    p_ASCAT_SM(i)=i/numel(tt__SM);
end
t_max=tt__SM(min(find(p_ASCAT_SM>0.95)))

t_avg=[]; k=1;
 for i=(mapping_start_doy-start_doy)*48+1:(mapping_end_doy-start_doy+1)*48
     t=(m_SM(a1:a2,b1:b2,i));
     t_avg(k)=nanmean(t(:));
     k=k+1;
 end
% SM TS avg
figure(1)
filename = [ofp,study_area,'_TS_SM_all_',num2str(no_data_mask(1)),num2str(no_data_mask(2)),num2str(no_data_mask(3)),...
    num2str(no_data_mask(4)),num2str(no_data_mask(5)),'_city.gif'];
k=1;
for i=(mapping_start_doy-start_doy)*48+1:run_hour*2:(mapping_end_doy-start_doy+1)*48
    if k==1
        
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

    else
        hold on
        plot(title_num((mapping_start_doy-start_doy)*48+1:(mapping_end_doy-start_doy+1)*48),t_avg, 'linewidth',1, 'color','k')
        ylim([0.2 0.35])
        datetick('x','mm-dd')
        set(gca,'fontweight','bold', 'fontsize',20)
        ylabel('Soil Moisture (m^3/m^3)','fontsize',30)
        ay=get(gca,'YtickLabel');
        set(gca, 'YtickLabel', ay,'fontsize',25)
        grid on       

        set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
        title({'Khartoum (Urban+Agriculture)','-Soil Moisture'},'fontsize',35)
        %title('Khartoum (Urban)-Soil Moisture','fontsize',35)

        if t_max<=t_avg(k*run_hour*2)
            plot(title_num(i),t_avg(k*run_hour*2),'ko','markerfacecolor',[0.8 0.3 0.3],'markersize',10)
        else
            plot(title_num(i),t_avg(k*run_hour*2),'ko','markerfacecolor',[0.3 0.3 0.8])
        end
    end
    box on
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    
    k=k+1;
end
