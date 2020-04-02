function all_metrics=bjerklie_suiteUNC(Q_act,finalAMHG)

   


%metrics
Qact=Q_act;% (:,columnindex==columncheck);
Q_reachavg=finalAMHG;
resid=((Qact-Q_reachavg)./Qact)*100;
sq=resid.^2;
SSE=sum(sq);
meansse=SSE/length(sq);
rrmse=sqrt(meansse);
relresid=resid;
mean_relresid=mean(relresid);
std_relresid=std(relresid);


resid=((Qact-Q_reachavg));
logresid=(log10(Qact)-log10(Q_reachavg));
sq=resid.^2;
SSE=sum(sq);
meansse=SSE/length(sq);
rmse=sqrt(meansse);
resids=resid;
mean_resid=mean(resids);
std_resid=std(resids);
mean_logresid=mean(logresid);
std_logresid=std(logresid);


% resid_mean=(Q_reachavg)-mean(Qact);
% resid=((Qact-Q_reachavg));
% sq_mean=resid_mean.^2;
% SSE_mean=sum(sq_mean);
% sq=resid.^2;
% SSE=sum(sq);
% error1=SSE_mean/SSE;
% m=2;%three fitting parameters-a/b, slope, int
% error2=(2*m)/size(Qact,2);
% MSC=log(error1-error2);


all_metrics=[rmse,rrmse,mean_resid,std_resid,mean_logresid,...
    std_logresid,mean_relresid,std_relresid];%num_time_removed,percent_cross_remain,


end
