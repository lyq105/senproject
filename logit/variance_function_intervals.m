function [m_inter,s_inter,lp_inter,varm,vars,sigma] = variance_function_intervals(x,n,r,nnum,datatype,cf_level,p)
%  计算fisher信息矩阵 
%   数据（x，n，r） theta=（theta1，theta2）
%   输入：
%             x :  刺激水平
%             n :  相同刺激水平试验次数
%             r :  相同刺激水平响应次数
%          nnum :  刺激水平个数 
%         theta :  待估计参数
%      datatype :  'norm'代表正态分布，'logistic'代表logit分布
%      cf_level :  置信度
%             p :  发火点概率
%   输出：
%       m_inter :  参数m的置信区间
%       s_inter :  参数s的置信区间
%      lp_inter :  Lp的置信区间
%          varm :  mu的方差
%          vars :  sigma的方差
%         sigma :  sigma的无偏估计
[theta_e,fval]= maximum_likelihood_estimates( x,n,r,nnum,datatype)
   
mu = theta_e(1);
sigma_e = theta_e(2)*pi/sqrt(3);

%%%!!! 这里存在问题 Analysis of Sensitivity Tests 文章的公式(6)与软件中的计算方法不同。
%beta = 1.6 /n^(2/3)     %% 文章中的取法
beta = 2/3.0               %% 软件中样本量为8时的取法

sigma = (1+beta)*sigma_e;

varm = 2.5*sigma^2 /nnum;
vars = 3.2*sigma^2 /nnum;
%vars = 1.2*sigma^2/n^(2/3)

cf = 1 - 0.5*(1 - cf_level);
cm1 = mu - sqrt(varm)*icdf(datatype,cf,0,sqrt(3)/pi);
cm2 = mu + sqrt(varm)*icdf(datatype,cf,0,sqrt(3)/pi);
cs1 = sigma - sqrt(vars)*icdf(datatype,cf,0,sqrt(3)/pi);
cs2 = sigma + sqrt(vars)*icdf(datatype,cf,0,sqrt(3)/pi);
m_inter = [cm1',cm2'];
s_inter = [cs1',cs2'];

lp_inter = [];
for i=1:length(p)
	zp = icdf(datatype,p(i),0,sqrt(3)/pi);
	lp1 =  mu + zp*sigma - sqrt(varm + zp*zp *vars)*icdf(datatype,cf,0,sqrt(3)/pi);
	lp2 =  mu + zp*sigma + sqrt(varm + zp*zp *vars)*icdf(datatype,cf,0,sqrt(3)/pi);
	lp_inter = [lp_inter,lp1',lp2'];
end

end

