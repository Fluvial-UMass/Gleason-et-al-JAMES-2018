
clear all
close all
clc

%Code to read this Nile_Landsat_RivWidth matrix

% Description:
% This matrix (Nile_Landsat_RivWidth.mat) documents cross-sectional flow widths
% along the Nile section ranging from north Faiyum to Samalut, Egypt,
% on 90 Landsat snapshots acquired during 1984-2015.

% River inundation extents were mapped by a local adaptive water extraction algorithm
% (please refer to Li and Sheng 2012, Wang et al. 2014 for details);

% River flow widths were then calculated by RivWidth version4 (please refer
% to Pavelsky et al. (2008) for details).

% The matrix is in a dimension of:
% 5272 rows (cross sections) x 7 columns (attributes) x 90 bands (dates)

% The 5272 cross sections are ordered longitudinally from north to south

% The 7 attributes include (from left to right):
% 1. Cross-section ID (1 to 5272)
% 2. Longitudinal distance (in meters) from the start (northern-most) cross section
% 3. X coordinate values (Easting in UTM Zone 36, in meters)
% 4. Y coordinate values (Northing in UTM Zone 36, in meters)
% 5. Cross section's reach id
% 6. Channel number at this cross section (a number > 1 indicates braided channels at this cross section)
% 7. Flow width (flow width at this cross section in meters)
% 8. Aquisition year
% 9. Acquisition month
% 10. Acquisition day

% poolobj=gcp('nocreate');
% if isempty(poolobj)
% parpool;
% end



% SAMPLE CODE TO MANIPULATE THIS MATRIX:----
load 'D:\Research\Projects\Nile\Nile for Jida\Nile\Nile_Landsat_RivWidth_mar17_2016.mat';
load Nile_Q_simulations.mat
gauge_in=xlsread('D:\Research\Projects\Nile\Nile for Jida\Nile\observed nile 1973 to 1984.xlsx');

t1=clock;
runit=1;
cross_sectional_var_filter=1;%remove low b cross sections
write_output=1;
plot_output=1;
meanchoice=4; %1 for model, 2 for observed, 3 for P/ET, 4 for ds/dt
n_cross_touse=15;
%----------------------------------

%------------------------
%set low and high filters for maximum and minimum reaosnable discharge.
%This is critical in eliminating false convergence of the GA outside of
%reaosnable ranges
lowfilter=100;%(10^Qmean_lo)/filterfactor;%cms
hifilter=20000;%(10^Qmean_hi)*filterfactor;%cms
%----------------

%set GA parameters.
blo=.05;
bhi=.7;

qweight=1;
wweight=1;
pop=20;
gen=5000;
%----------------

%set parameters for Qc-------------
if meanchoice==1
    Qmean_wet=(12074);%cms
    Qmean_dry=(10684);%cms
    meanstring='model mean';
elseif meanchoice==2
    Qmean_wet=(1440);%cms
    Qmean_dry=(1120);%cms
    meanstring='gauge mean';
elseif meanchoice==3
    Qmean_wet=(6000);%cms
    Qmean_dry=(2000);%cms
    meanstring='P ET mean';
elseif meanchoice==4
    Qmean=3666;%cms
    meanstring='dsdt mean';
end
%assuming a 70% variation  in mean Q
qchi=log10(Qmean+(Qmean*0.7));
qclo=log10(Qmean-(Qmean*0.7));
%------------------------------------

% name the output
string1='D:\Research\Projects\Nile\Nile for Jida\Nile\textfiles\Nile_amhg_dsdt_oct9_changes_multiobj2';
string2=[num2str(n_cross_touse), 'cross_'];
string3=[num2str(gen), 'gen_'];
string4='';
string5=meanstring;
outfile1=[string1,string2,string3,string4,string5,'amhg1.txt'];
outfile2=[string1,string2,string3,string4,string5,'amhg2.txt'];
outfile3=[string1,string2,string3,string4,string5,'amhg3.txt'];
outfile4=[string1,string2,string3,string4,string5,'amhg4.txt'];
outfile_date=[string1,string2,string3,string4,string5, '_date.txt'];
outfile_rmse=[string1,string2,string3,string4,string5, '_widthrmse.txt'];
        

%prepare input data
W=reshape(Nile_Landsat_RivWidth(:,7,:),size(Nile_Landsat_RivWidth,1),size(Nile_Landsat_RivWidth,3));
dates_input_year=squeeze(Nile_Landsat_RivWidth(:,8,:));
dates_year=dates_input_year(1,:)';

dates_input_month=squeeze(Nile_Landsat_RivWidth(:,9,:));
dates_month=dates_input_month(1,:)';

dates_input_day=squeeze(Nile_Landsat_RivWidth(:,10,:));
dates_day=dates_input_day(1,:)';

dates=horzcat(dates_year,dates_month,dates_day);



%rmeoving  dates with W<92m
index_of_removal=[1:7,8:37,39:72,74:75,77:78,80:91];
W=W(:,index_of_removal);
dates_year=dates_year(index_of_removal);
dates_month=dates_month(index_of_removal);
dates_day=dates_day(index_of_removal);
dates=dates(index_of_removal,:);

[dates_inorder,dateindex]=sort(datenum(double(dates)));
W=W(:,dateindex);  %           
dates_month_inorder=dates_month(dateindex);
all_dates=dates_inorder;

%----------------

%now, break them up into predefiend reaches based on inflows/outflows
reachindex_in=squeeze(Nile_Landsat_RivWidth(:,5,:));
reachindex=reachindex_in(:,1);
%--------------------

%now, analyzie width variability of the data
[W,width_nospace,flag,filter_index]=remove_low_width_crossections(W,cross_sectional_var_filter);
reachindex=reachindex(filter_index,:);
reach_ids=unique(reachindex(reachindex>=0));
%looks good!

      %  

%----------------
if runit==1
    
    %for every reach, we need to
    %     1) subsample the cross sections therein to ease computational burden
    %     2) set AMHG parameters by season
    %     3) run the meothod by season and then combine the hydrographs
    %     final_amhg_wet=NaN(60,reachindex(1));d
    %     final_amhg_dry=NaN(60,reachindex(1));
    %     finalrmse_wet=NaN(n_cross_touse,reachindex(1));
    %     finalrmse_dry=NaN(n_cross_touse,reachindex(1));
    %     all_amhg=zeros(87,reachindex(1));
    %     alldates_out=zeros(87,reachindex(1));
    
    all_amhg1=NaN(length(dates_inorder),reachindex(1));
    %[JW]: there are 28 unique reaches (0 to 27), but this demension is
    %[JW]: only 87 dates x 27 reaches. 
    all_amhg2=NaN(length(dates_inorder),reachindex(1));
    all_amhg3=NaN(length(dates_inorder),reachindex(1));
    all_amhg4=NaN(length(dates_inorder),reachindex(1));
    alldates_out=NaN(length(dates_inorder),reachindex(1));
    
    for reachnum=0:reachindex(1)
        
        disp(['reach number ', num2str(reachnum)])
       
        
        reachW=W(reachindex==reachnum,:);
        
        
        %set number of cross sections for a test or subset. Set this to
        %size(width,1) to use all cross sectins, but be prepared to wait many
        %hours in that case
                
         if size(reachW>15) %[JW]: Colin, is this to test the cross section number?
            %[JW]: if so, should it be size(reachW, 1)>15
            randindex=randi(size(reachW,1),n_cross_touse,1);
            sampled_W=reachW(randindex,:);
            
        elseif size(reachW)<5 %[JW]: Similarly, should this be size(reachW, 1)<15
            continue
        else
            sampled_W=reachW;
        end
        
        %rmeove low width cross sections
        
        [sampled_W,junk,flag1,junk]=remove_low_width_crossections(sampled_W,cross_sectional_var_filter);
        
        disp('     ')
        if flag1 ==1
            disp(['too low b!'])
            continue
        end
        
        pint=calculate_pint(sampled_W);
        
        if pint<0.15
            flag1=1;
            disp(['reach ', num2str(reachnum), ' has pint = ', num2str(pint), ' and thus no AMHG'])
        end
        
        if flag1 ==1
            disp(['no AMHG here'])
            continue
        end
        
        
        
        %dlmwrite(['D:\Nile JAMES\textfiles\reach no',num2str(reachnum),'.txt'],sampled_W)
        
        
        %set AMHG parameters, following Gleason and Wang 2015's defintion of
        %the slope and interecept of AMHG. Slope = mode(1/logQmean),
        %int=mode(logWmean)/mode(logQmean). The wider the range in slope adn
        %intercept, the more variable the pivot poitn of rating curves and the
        %more the need for longer or multiple GAs  
        
          
        abslopelo=(-1/(qclo));
        abslopehi=(-1/(qchi));
        
        wclo=mean(mode(round(log10(reachW)*10)/10))*0.99;
        wchi=mean(mode(round(log10(reachW)*10)/10))*1.01;
        
        abintlo=-wclo*abslopehi;%calculating slope from int
        abinthi=-wclo*abslopelo;
        
        avect=[qclo:(qchi-qclo)/50:qchi]';%make a 50 unit search space
        bvect=[wclo:(wchi-wclo)/50:wchi]';
        qcwcvect=horzcat(avect,bvect);
        
        
        %----------------
        
        
 
        
        
        if flag1~=1
            
            t3=clock;
            [meanQout,medianQout]=invoke_AMHG_ngsa(size(sampled_W,1),sampled_W,...
                lowfilter,hifilter,bhi,blo,qclo,qchi,wclo,wchi,pop,gen,qcwcvect,qweight,wweight);
            
            
            t4=clock;
            elapsed=datevec(datenum(t4)-datenum(t3));
            disp(['reach number ', num2str(reachnum),' took ', num2str(elapsed(4))...
                ' hours, ', num2str(elapsed(5)), ' min, and ', num2str(elapsed(6)), ' seconds' ])
            disp('    ')
            
        end
        
        if flag1==1
            Q1=NaN(1,size(W,2));
            Q2=NaN(1,size(W,2));
            Q3=NaN(1,size(W,2));
            Q4=NaN(1,size(W,2));
            
        else
            Q1=nanmean(meanQout);
            Q2=nanmedian(meanQout);
            Q3=nanmean(medianQout);
            Q4=nanmedian(medianQout);
        end
        
        
        %stitch them together
        all_amhg1(1:length(Q1),reachnum+1)=(Q1');  
        all_amhg2(1:length(Q2),reachnum+1)=(Q2');
        all_amhg3(1:length(Q3),reachnum+1)=(Q3');
        all_amhg4(1:length(Q4),reachnum+1)=(Q4');
        
        alldates_out(1:length(Q1),reachnum+1)=all_dates;
        
        
        
        if write_output==1
            dlmwrite(outfile1,all_amhg1)
            dlmwrite(outfile2,all_amhg2)
            dlmwrite(outfile3,all_amhg3)
            dlmwrite(outfile4,all_amhg4)
            dlmwrite(outfile_date,alldates_out,'precision',8)
            
        end
        
        
        
        
    end
    
    
    
    
    t2=clock;
    elapsed=datevec(datenum(t2)-datenum(t1));
    disp(['total time ' num2str(elapsed(4))...
        ' hours, ', num2str(elapsed(5)), ' min, and ', num2str(elapsed(6)), ' seconds' ])
    
    
end



if plot_output==1
    figure;
    allamhgin1=dlmread(outfile1);
    allamhgin2=dlmread(outfile2);
    allamhgin3=dlmread(outfile3);
    allamhgin4=dlmread(outfile4);
    alldatesin=dlmread(outfile_date);
    
    monthly_gauge_in=xlsread('D:\nile\observed monthly nile 1973 to 1984.xlsx');
    gauge_dates=datenum([monthly_gauge_in(:,1:2),ones(length(monthly_gauge_in),1)]);
    
    plot(gauge_dates,monthly_gauge_in(:,4),'k', 'LineWidth',5)
    hold on
    plot(alldatesin,allamhgin1,'rx')
    
    plot(alldatesin,allamhgin2,'bx')
    
    plot(alldatesin,allamhgin3,'gx')
    
    plot(alldatesin,allamhgin4,'cx')
    
    datetick('x',2,'keepticks')
    xlabel('Date','FontSize',16)
    ylabel('Q cms','FontSize',16)
    legend('Gauge monthly','AMHG mean mean','AMHG median mean','AMHG mean median','AMHG median median',...
        'Location','NorthWest')
    
end

nanmedian(nanmedian(allamhgin4))
nanmean(nanmean(allamhgin4))

nanmedian(nanmedian(allamhgin3))
nanmean(nanmean(allamhgin3))

plot(allamhgin4)

% alldates_out=repmat(all_dates,1,28);
%  dlmwrite(outfile_date,alldates_out,'precision',8)

