function [solutions1,flag]=AMHG_initilization_ratio_UMass(qclo,qchi,wc,width1,width2,lowfilter,hifilter,bhi,blo)


flag=0;

% 
% wc=mean([wchi,wclo]);
% width1=sampled_W(1,:);
% width2=sampled_W(2,:);

%below is a solution for b wrt AMHG. the search space given these
%parameters is defined by these b values, which solve for the minimum b
%needed to make output Q pass Q filters given input widths. the first b in
%each case solves for the hi fitler, the second for the low filter.

% width1=log10(width1);
% width2=log10(width2);


solutions1=NaN(100,2);
bvect=blo:(bhi-blo)/50:bhi;
qcvect=qclo:(qchi-qclo)/50:qchi;


% for i=1:100
%     current_b=bvect(i);
%     for j=1:100
%         current_qc=qcvect(j);
%         plot(   (    log10(width1)+(current_b.*current_qc)-wc     )   ./(current_b),'b')
%         hold on
%         
%     end
% end
% 
% plot(1:90,ones(90,1).*log10(lowfilter))
% plot(1:90,ones(90,1).*log10(hifilter))
% ylim([0 5])


% 

count1=0;
count2=0;


for j=1:50
    current_qc=qcvect(j);
    
    for i=1:50
        current_b=bvect(i);
        for k=1:50
            
         current_b2=bvect(k);
         
        hivect = (    log10(width1)+(current_b.*current_qc)-wc     )   ./(current_b);
        hivect2 = (    log10(width2)+(current_b2.*current_qc)-wc     )   ./(current_b2);
        if all(hivect>log10(lowfilter)) && all(hivect<log10(hifilter))&& all(hivect2<log10(hifilter))&& all(hivect2>log10(lowfilter))
            count1=count1+1;
            solutions1(count1,1:3)=[current_b,current_b2,current_qc];
               
        end
        
        end
    end
end



if all(isnan(solutions1))
    flag=1;
end




