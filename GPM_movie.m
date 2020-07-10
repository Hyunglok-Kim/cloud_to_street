clc; clear

addpath('/sfs/qumulo/qproject/hydrosense/shape_files/sudan/')
addpath('/sfs/qumulo/qproject/hydrosense/shape_files/moza/')
shp_file_dir='/sfs/qumulo/qproject/hydrosense/shape_files/';
ifp=['/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_01/flood_movie/'];
ofp='/project/hydrosense/to_OUTSIDE/';

year_=2019;
study_area='sudan'
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
    start_doy=210; end_doy=250; % expected dates of flooding +/-
    mapping_start_doy=220; mapping_end_doy=244;
    
    % h-box scale
    a1=20;a2=55;
    b1=15; b2=45;
    r1=1.2; r2=2; title_position=0.92;
    
    %Khartoum City Flooding: (focused on urban flooding)
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

hh_start_doy=(start_doy-1)*48+1;
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
% load data
clc
%t_P=load([iifp,'GPM_PR_',num2str(year_),'_',num2str(start_doy),'_',num2str(end_doy),'.mat']);
t_P=load([iifp,'GPM_PR_',num2str(year_),'.mat']);
GPM_P=t_P.GPM_PR;
GPM_P=GPM_P(:,:,hh_start_doy:hh_end_doy);
%%
% check spatail averaged P data
 t_avg=[]; k=1;
 for i=1:size(GPM_P,3)
     t=(GPM_P(a1:a2,b1:b2,i));
     t_avg(k)=nanmean(t(:));
     k=k+1;
 end

plot(title_num,t_avg,'-')
axis ij
datetick('x','yy-mm-dd')
ylabel('Precipitation (mm/hr)','fontsize',15)
grid on
set(gca,'fontweight','bold', 'fontsize',15)
%ylim([0.15 0.35])
title('Khartoum 2019 - Hyetograph','fontsize',20)

set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
myStyle = hgexport('factorystyle');
myStyle.Format = 'tiff';
%myStyle.Width = 1;
%myStyle.Height = 1;
myStyle.Resolution = 100;
myStyle.Units = 'inch';
myStyle.FixedFontSize = 20;

%hgexport(gcf,['/project/hydrosense/to_OUTSIDE/P.tiff'],myStyle,'Format','tiff')
close(gcf)
%%
% calculate 95% of P
a1=24;a2=27;
b1=25;b2=28;

%a1=23;a2=45;
%b1=20;b2=38;

t_GPM=GPM_P(a1:a2, b1:b2,:);
tt_GPM=[];
for i=1:size(t_GPM,3)
    t=t_GPM(:,:,i);
    tt_GPM(i)=nanmean(t(:));
end
tt_GPM=sort(tt_GPM);
p_GPM=[];
for i=1:numel(tt_GPM)
    p_GPM(i)=i/numel(tt_GPM);
end
t_max=tt_GPM(min(find(p_GPM>0.95)));

t_avg=[]; k=1;
 for i=(mapping_start_doy-start_doy)*48+1:run_hour*2:(mapping_end_doy-start_doy+1)*48
     t=(GPM_P(a1:a2,b1:b2,i));
     t_avg(k)=nanmean(t(:));
     k=k+1;
 end
 plot(title_num((mapping_start_doy-start_doy)*48+1:run_hour*2:(mapping_end_doy-start_doy+1)*48),t_avg, 'linewidth',1, 'color',[0.1 0.5 0.8])
 axis ij
 datetick('x','mm-dd')
 set(gca,'fontweight','bold', 'fontsize',20)
 ylabel('Precipitation (mm/hr)','fontsize',30)
 grid on
 ylim([0 2.5])
 ay=get(gca,'YtickLabel');
 set(gca, 'YtickLabel',ay,'fontsize',25)
 
 %title({'Khartoum (Urban+Agriculture)','-Precipitation'},'fontsize',35)
 title('Khartoum (Urban)-Precipitation','fontsize',35)
 set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
