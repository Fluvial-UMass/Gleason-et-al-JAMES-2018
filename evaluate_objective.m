function f = evaluate_objective(x, M, V,width1,width2,qweight,wweight,lowfilter,hifilter);

% function f = evaluate_objective(x, M, V)
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables.
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input. Make sure that the
%
% An example objective function is given below. It has two six decision
% variables are two objective functions.

% f = [];
% %% Objective function one
% % Decision variables are used to form the objective function.
% f(1) = 1 - exp(-4*x(1))*(sin(6*pi*x(1)))^6;
% sum = 0;
% for i = 2 : 6
%     sum = sum + x(i)/4;
% end
% %% Intermediate function
% g_x = 1 + 9*(sum)^(0.25);
%
% %% Objective function two
% f(2) = g_x*(1 - ((f(1))/(g_x))^2);

% Kursawe proposed by Frank Kursawe.
% Take a look at the following reference
% A variant of evolution strategies for vector optimization.
% In H. P. Schwefel and R. Männer, editors, Parallel Problem Solving from
% Nature. 1st Workshop, PPSN I, volume 496 of Lecture Notes in Computer
% Science, pages 193-197, Berlin, Germany, oct 1991. Springer-Verlag.
%
% Number of objective is two, while it can have arbirtarly many decision
% variables within the range -5 and 5. Common number of variables is 3.
f = [1,1];
% Objective function one


b1=x(1);
b2=x(2);
logqc=(x(3));
logwc=(x(4));


logq1=(log10(width1)+(b1.*logqc)-logwc)./b1;
logq2=(log10(width2)+(b2.*logqc)-logwc)./b2;

width2_invert=(x(2).*logq1)-(x(2)*logqc)+logwc;
width2_invert=10.^width2_invert;



Qresid=((10.^logq1-10.^logq2));%./(logq1+logq2));
Qsq=Qresid.^2;
QSSE=sum(Qsq);
Qmeansse=QSSE/length(Qsq);
QRMSE=sqrt(Qmeansse);




Wresid=((width2-width2_invert));  
Wsq=Wresid.^2;
WSSE=sum(Wsq);
Wmeansse=WSSE/length(Wsq);
WRMSE=sqrt(Wmeansse);

% Decision variables are used to form the objective function.
f(1) = QRMSE;%qweight*
% Decision variables are used to form the objective function.
f(2) = WRMSE;%wweight*

if any(10.^logq1<lowfilter) || any(10.^logq2<lowfilter)|| any(10.^logq1>hifilter)|| any(10.^logq2>hifilter)
    f(1)=1e15;
    f(2)=1e15;
end


% Check for error
if length(f) ~= M
    error('The number of decision variables does not match you previous input. Kindly check your objective function');
end