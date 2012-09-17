function [theta_e,fval]=maximum_likelihood_estimates( x,n,r,nnum,datatype)
% 计算极大似然估计
%   数据（x，n，r） theta=（theta1，theta2）
%   x :  刺激水平
%   n :  相同刺激水平试验次数
%   r :  相同刺激水平响应次数
%   nnum : 刺激水平个数 
%   theta：待估计参数
%   datatype: 感度数据类型 'norm':正态分布，'logistic'：Logistic分布

[theta_e,fval] = fminsearch(@(theta)-log(likelihood_function(x,n,r,nnum,theta,datatype)),[mean(x),var(x)]);
fval = - fval;
end 
