function [meanQout,medianQout]=invoke_AMHG_ngsa(numcrosssections,input_width,...
    lowfilter,hifilter,bhi,blo,qclo,qchi,wclo,wchi,pop,gen,qcwcvect,qweight,wweight)

% % lowfilter=lowfilter_dry;
% % hifilter=hifilter_dry;
% % input_width=sampled_W_dry;
% % width1=sampled_W_dry(1,:);
% % width2=sampled_W_dry(2,:);
% % qcwcvect=qcwcvect_dry;
% % numcrosssections=size(sampled_W_dry,1);

%run the GA----------------------------------------------------
%run the GA----------------------------------------------------

combinations=nchoosek(1:numcrosssections,2);
meanQout=NaN(length(combinations),size(input_width,2));
medianQout=meanQout;
misscount=0;
wc=mean(qcwcvect(:,2));

parfor p=1:length(combinations)
    
    k=combinations(p,1);
    i=combinations(p,2);
    width1=input_width(k,:);
    width2=input_width(i,:);
    %
    %check for soultions based on limits here
    % [blox1,blox2,flag]=AMHG_initilization_ratio_UNC(qcwcvect,width1,width2,lowfilter,hifilter,bhi,blo);
    [solutions,flag]= AMHG_initilization_ratio_UMass(qclo,qchi,...
        wc,width1,width2,lowfilter,hifilter,bhi,blo);
    
    
    if flag==1
        %there aren't any solutions within limits
        %this is why holidt is NaN
        misscount=misscount+1;
        continue
    end
    t1=datenum(clock);
    paretosolutions=nsga_2(pop,gen,width1,width2,solutions,qcwcvect,bhi,qweight,wweight,lowfilter,hifilter);
    holdit=NaN(pop,size(input_width,2));
    elapsed=datevec(datenum(clock)-datenum(t1))
    
    paretosolutions_rank1=paretosolutions(paretosolutions(:,7)==1,:);
    paretosolutions_rank1=paretosolutions_rank1(paretosolutions_rank1(:,6)<90,:) %3 pixel width error
    paretosolutions_rank1=paretosolutions_rank1(paretosolutions_rank1(:,5)<500,:) %500m3/s flow error
    
    if isempty(paretosolutions_rank1)
        continue
    end
    
    for j=1:size(paretosolutions_rank1,1)
        b1=paretosolutions(j,1);
        b2=paretosolutions(j,2);
        logqc=(paretosolutions(j,3));
        logwc=(paretosolutions(j,4));
        
        q1=10.^((log10(width1)+b1.*logqc-logwc)./b1);
        q2=10.^((log10(width2)+b2.*logqc-logwc)./b2);
        qout=mean([q1;q2]);
        
        
        holdit(j,:)=qout;
    end
    
    
    meanQout(p,:)=nanmean(holdit);
    medianQout(p,:)=nanmedian(holdit);
end
%     figure;
%         plot(nanmean(holdit))
%         hold on;
%         plot(nanmedian(holdit))
%         plot(holdit')
end




%-----------------------------------------------------------------------


%  -----------------------------------------------------------------------------