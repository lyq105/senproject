function [m_inter,s_inter,lp_inter] = asymptotic_confidence_intervals(x,n,r,nnum,datatype,cf_level,p)
%  计算fisher信息矩阵 
%   数据（x，n，r） theta=（theta1，theta2）
%   输入：
%             x :  刺激水平
%             n :  相同刺激水平试验次数
%             r :  相同刺激水平响应次数
%          nnum :  刺激水平个数 
%         theta :  待估计参数
%      datatype :  感度数据类型 0:正态分布，1：Logistic分布
%      cf_level :  置信度
%             p :  发火点概率
%   输出：
%       m_inter :  参数m的置信区间
%       s_inter :  参数s的置信区间
%      lp_inter :  Lp的置信区间
%
[theta_e,fval]= maximum_likelihood_estimates( x,n,r,nnum,datatype)
   
[i00,i01,i11] = fisher_information_matrix(x,n,r,nnum,theta_e,datatype);

v = 0.5*atan(2*i01/(i00-i11));

d_cf = chi2inv(cf_level,1);
a = sqrt(2 * d_cf / (i00 + i11 + sqrt((i00 - i11)^2 + (2 * i01)^2)));
b = sqrt(2 * d_cf / (i00 + i11 - sqrt((i00 - i11)^2 + (2 * i01)^2)));

alpha = 0:1e-5:2*pi;

m_inter =[];
s_inter =[];
lp_inter=[];

for i=1:length(d_cf)
	m = theta_e(1) + a(i)*cos(alpha)*cos(v) - b(i)*sin(alpha)*sin(v);
	s = theta_e(2) + a(i)*cos(alpha)*sin(v) + b(i)*sin(alpha)*cos(v);
		
	m_inter = [m_inter;[min(m),max(m)]];
	s_inter = [s_inter;[min(s),max(s)]];

	lp_in_t = [];
	for j=1:length(p)
		lp =  m + icdf(datatype,p(j),0,1)*s;
		lp_in_t = [lp_in_t,min(lp),max(lp)];
	end
	lp_inter = [lp_inter;lp_in_t];
end

end

