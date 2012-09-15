function ret=likelihood_function( x,n,r,nnum,theta,datatype)
%   两种感度数据的似然函数
%   数据（x，n，r） theta=（theta1，theta2）
%   x :  刺激水平
%   n :  相同刺激水平试验次数
%   r :  相同刺激水平响应次数
%   nnum : 刺激水平个数 
%   theta：待估计参数
%   datatype: 感度数据类型 0:正态分布，1：Logistic分布

ret = 1;
for i =1:nnum
      z = (x(i) - theta(1)) / theta(2);
      ret = ret* combinatorialnumber(n(i),r(i)) * (cdf(datatype,z,0,1)^(r(i))) * ((1-cdf(datatype,z,0,1))^(n(i) - r(i)));  
end

end

