clc
clear all
close all




% modellat= ncread('D:/Nile/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile_amhg.nc4','lat');
% modeltime=ncread('D:/Nile/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile_amhg.nc4','time');
% modelQ=   ncread('D:/Nile/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile_amhg.nc4','discharge');
% modelQ(modelQ<50)=NaN;
% 
modellat= ncread('D:/Nile JAMES/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile.nc4','lat');
modeltime=ncread('D:/Nile JAMES/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile.nc4','time');
modelQ=   ncread('D:/Nile JAMES/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile.nc4','discharge');
modelQ_monthly=ncread('D:/Nile JAMES/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile_month.nc4','discharge');
modelQ_monthly_time=ncread('D:/Nile JAMES/Global-WFDEI-Extended-historical-varsoc-discharge30min_nile_month.nc4','time');
modelQ(modelQ<50)=NaN;
modelQ_monthly(modelQ_monthly<50)=NaN;

%pull flow only from our reach, which covers latitude 28.25 to 29.75
lat_index=(modellat>= 28.25 & modellat<= 29.75) ;
meanQ=squeeze(nanmean(nanmean(modelQ(:,lat_index,:))));
offset=((datenum('Jan 1 1901')-datenum('January 1 0000')));
modeltimeout=datenum(modeltime+offset);
[sortedtime,index]=sort(modeltimeout);
tuned_sortedmeanQ=meanQ(index);
tuned_sorteddate=sortedtime;

%monthly flow for our reach
mean_month_Q=squeeze(nanmean(nanmean(modelQ_monthly(:,lat_index,:))));
monthtime=datenum(modelQ_monthly_time+offset);
[sortedtime_month,index]=sort(monthtime);
mean_month_Q=mean_month_Q(index);
monthtime=sortedtime_month;



load Nile_Q_simulations.mat 
%pre amhg simulations

oldmodelQ=Nile_Q_simulations(:,6);
oldmodeldate=datenum(double(horzcat(Nile_Q_simulations(:,1),Nile_Q_simulations(:,2),Nile_Q_simulations(:,3))));


gauge1='D:\Nile\textfiles\Nile_amhg_dsdt_oct9_changes_multiobj215cross_5000gen_dsdt meanamhg1.txt';
gauge2='D:\Nile\textfiles\Nile_amhg_dsdt_oct9_changes_multiobj215cross_5000gen_dsdt meanamhg2.txt';
gauge3='D:\Nile\textfiles\Nile_amhg_dsdt_oct9_changes_multiobj215cross_5000gen_dsdt meanamhg3.txt';
gauge4='D:\Nile\textfiles\Nile_amhg_dsdt_oct9_changes_multiobj215cross_5000gen_dsdt meanamhg4.txt';

outfile_date='D:\Nile\textfiles\for yoshi dates.txt';

gauge_Discharge1=dlmread(gauge1);
gauge_Discharge2=dlmread(gauge2);
gauge_Discharge3=dlmread(gauge3);
gauge_Discharge4=dlmread(gauge4);

meanAMHG=mean(cat(3,gauge_Discharge1,gauge_Discharge2,gauge_Discharge3,gauge_Discharge4),3);


fid=fopen(outfile_date);
alldatesin=textscan(fid,'%s %s %s','Delimiter','.');
fclose(fid);
hmm1=alldatesin{1};
hmm2=alldatesin{2};
hmm3= alldatesin{3};

for i=1:length(alldatesin{1})
    dates_matlab(i,1)=datenum([hmm1{i},'/',hmm2{i},'/',hmm3{i}]);
end


monthly_gauge_in=xlsread('D:\nile\observed monthly nile 1973 to 1984.xlsx');
gauge_dates=datenum([monthly_gauge_in(:,1:2),ones(length(monthly_gauge_in),1)]);
gauge_flow=monthly_gauge_in(:,4);



% light blue, dark green, dark blue, light green, peach, light red
colorstring={[0 0 0.5],[0 1 0],[0 0 1],[ 0 0.5 0],[1 0.6 0.6],[0.8 0 0]};
symbolstring={'o','h','^','d','s','v'};
%---------------------------------------


% % summary plot
% figure('Position',[100,100,1200,900])
% 
% plot(gauge_dates,(gauge_flow),'k', 'LineWidth',3)
% hold on
% plot(oldmodeldate,(oldmodelQ),'Color',[0 0.5 0],'LineWidth',2)
% plot(tuned_sorteddate,(tuned_sortedmeanQ),'r','LineWidth',2)
% plot(dates_matlab,(meanAMHG),'bx')
% 
% datetick('x',2,'keepticks')
% xlabel('Date','FontSize',16)
% ylabel('Discharge (m^3/s)','FontSize',16)
% lgd=legend('Historical gauge','Model before AMHG tuning','Model after AMHG tuning','AMHG Q estimations','Location','NorthWest');
% lgd.FontSize=16;
% title('Tuning the PCR-GLOBWB model with AMHG over the lower Nile','FontSize',16)
% set(gca,'FontSize',16,'LineWidth',2)
% %---------------------------------------

% summary plot
figure('Position',[100,100,1200,900])

plot(gauge_dates,(gauge_flow),'k', 'LineWidth',3)
 hold on
 plot(oldmodeldate,(oldmodelQ),'Color',[0 0.5 0],'LineWidth',2)
 plot(dates_matlab,(meanAMHG),'bx')
 plot(tuned_sorteddate,(tuned_sortedmeanQ),'r','LineWidth',2)
xticks(min(gauge_dates):1095.75:max(tuned_sorteddate))
datetick('x',2,'keepticks')
xtickangle(45)
xlabel('Date','FontSize',16)
ylabel('Discharge (m^3/s)','FontSize',16)
lgd=legend('Historical gauge','Model before AMHG tuning','AMHG Q estimations','Model after AMHG tuning','Location','NorthWest');
lgd.FontSize=16;
%title('Tuning the PCR-GLOBWB model with AMHG over the lower Nile','FontSize',16)
set(gca,'FontSize',16,'LineWidth',2)
xlim([min(gauge_dates),max(tuned_sorteddate)])
ylim([0 3e4])
%---------------------------------------



%limit the historical gauge and tuned flow to the same number of months
date_start= monthtime(1);
date_end=monthtime(70);

comparison_model_Q=mean_month_Q(monthtime>=date_start & monthtime<date_end);
comparison_model_time=monthtime(monthtime>=date_start & monthtime<date_end);
comparison_gauge_Q=gauge_flow(gauge_dates>=date_start & gauge_dates<=date_end);
comparison_gauge_time=gauge_dates(gauge_dates>=date_start & gauge_dates<=date_end);
%---------------------------------------   


%plotting--------------------------------
figure('Position',[100,100,1200,900])
hold on;

plot(comparison_gauge_time,comparison_gauge_Q,[symbolstring{1},'-'],'Color','k','MarkerFaceColor','k');
plot(comparison_model_time,comparison_model_Q,[symbolstring{6},'-'],'Color',colorstring{6},'MarkerFaceColor',colorstring{6});

ylabel('Discharge (m^3/s)','FontSize',20)
xlabel('Date','FontSize',14)

xlim([min(comparison_model_time),max(comparison_model_time)])
xticks(min(comparison_model_time):180:max(comparison_model_time))
%plot(comparison_model_time-60,comparison_model_Q,['--'],'Color',colorstring{5},'MarkerFaceColor',colorstring{5},'LineWidth',1)
lgd=legend('Historical Gauge','Tuned Model','Location','NorthWest') %'Lagged Model',
lgd.FontSize=16;

%title('Tuned Model vs. Historical Gauge','FontSize',14)
datetick('x',12,'keepticks')
set(gca,'FontSize',16,'LineWidth',2)
%-----------------------------------------


%assessment metrics---------------------------
%[rmse,rrmse,mean_resid,std_resid,mean_logresid,std_logresid,mean_relresid,std_relresid]
allmetrics=bjerklie_suiteUNC(comparison_model_Q,comparison_gauge_Q)

%with a 2 month lag
comparison_model_Q=mean_month_Q(monthtime>=date_start & monthtime<date_end);
comparison_model_time=monthtime(monthtime>=date_start & monthtime<date_end);
comparison_gauge_Q=gauge_flow(gauge_dates>=date_start-60 & gauge_dates<=date_end-60);
comparison_gauge_time=gauge_dates(gauge_dates>=date_start-60 & gauge_dates<=date_end-60);
 
allmetrics_lagged=bjerklie_suiteUNC(comparison_model_Q(2:end-2),comparison_gauge_Q(2:end-2))
%-----------------------------------------



%plot AMHG against the model

date_start=datenum('Jan 1 2013');
date_end=datenum('Jan 1 2017');

comparison_model_Q=tuned_sortedmeanQ(tuned_sorteddate>=date_start & tuned_sorteddate<date_end);
comparison_model_time=tuned_sorteddate(tuned_sorteddate>=date_start & tuned_sorteddate<date_end);
comparison_AMHG_Q =meanAMHG(dates_matlab>=date_start &  dates_matlab<=date_end,:);
comparison_AMHG_time= dates_matlab( dates_matlab>=date_start &  dates_matlab<=date_end);



figure('Position',[100,100,1200,900])
hold on;

plot(comparison_model_time,comparison_model_Q,['-'],'Color',colorstring{6},'LineWidth',2);
plot(comparison_AMHG_time',nanmean(comparison_AMHG_Q,2),['-'],'Color',colorstring{1},'LineWidth',2);
plot(comparison_AMHG_time,comparison_AMHG_Q,['.'],'Color',colorstring{1});

% x1=comparison_AMHG_time;
% x=[x1;flipud(x1)];
% 
% y1=nanmean(comparison_AMHG_Q,2)+nanstd(comparison_AMHG_Q,0,2);
% y2=flipud(nanmean(comparison_AMHG_Q,2)-nanstd(comparison_AMHG_Q,0,2));
% y=[y1;y2];
% 
% fill(x,y,colorstring{1},'FaceAlpha',.2,'linestyle','none');

% errorbar(comparison_AMHG_time',nanmean(comparison_AMHG_Q,2),nanstd(comparison_AMHG_Q,0,2),'linestyle','none')

lgd=legend('Tuned model', 'AMHG','Location','NorthWest')
ylabel('Discharge (m^3/s)','FontSize',14)
xlabel('Date','FontSize',14)
% title('Tuned Model vs. Historical Gauge','FontSize',14)
datetick('x',12)
lgd.FontSize=19;
set(gca,'FontSize',16,'LineWidth',2)

%calculate sd for flow

nanmean(std(meanAMHG))/nanmean(nanmean(meanAMHG))






