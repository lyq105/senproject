function [varm,vars,covms] = cramer_rao_bounds( x,n,r,nnum,theta,datatype)
%  计算fisher信息矩阵 
%   数据（x，n，r） theta=（theta1，theta2）
%   输入：
%       x :  刺激水平
%       n :  相同刺激水平试验次数
%       r :  相同刺激水平响应次数
%       nnum : 刺激水平个数 
%       theta：待估计参数
%       datatype: 感度数据类型 0:正态分布，1：Logistic分布
%   输出：
%         Estimated Lower Bound of Variance of Mu，     varm
%         Estimated Lower Bound of Variance of Sigma    vars
%         Estimated CoVariance of Mu and Sigma          covms
%

[i00,i01,i11] = fisher_information_matrix(x,n,r,nnum,theta,datatype);

deti  =  i00*i11- i01*i01;
varm  =  i11/deti;	
vars  =  i00/deti*pi*pi/3;	
covms = -i01/deti*pi/sqrt(3);

end

