function [i00,i01,i11] = fisher_information_matrix( x,n,r,nnum,theta,datatype)
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
%       信息矩阵I = [i00,i01;i01,i11]
%

i00 = 0;
i01 = 0;
i11 = 0;
for i =1:nnum
      z = (x(i) - theta(1)) / theta(2);
			denominator = pdf(datatype,z,0,1)^2 / (cdf(datatype,z,0,1) * (1- cdf(datatype,z,0,1)) * theta(2)^2);
			i00 = i00 + n(i) * denominator;
			i01 = i01 + n(i) * z * denominator;
			i11 = i11 + n(i) * z^2 *denominator;
end

end

