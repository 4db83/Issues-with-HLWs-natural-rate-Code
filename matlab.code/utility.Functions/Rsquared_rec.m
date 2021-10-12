function [R2_exp R2_rec] = Rsquared_rec(u,y,Recession_indicator)
% Computes a Recession and Expension R-squared.
% db 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
US_rec 	= Recession_indicator;
y_			= demean(y);
% compute the 1-sum(I.*uhat.^2)/sum(I.*(r-rbar).^2) ratio
R2_exp	= 1 - mean(u(US_rec==0).^2)/mean(y_(US_rec==0).^2);
R2_rec	= 1 - mean(u(US_rec==1).^2)/mean(y_(US_rec==1).^2);

