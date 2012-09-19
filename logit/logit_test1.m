%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     logistic分布数据计算方法验证
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 数据
datatype = 'logistic';

x = [3.00000,6.07538,7.41994,9.00000,10.9945,12.0000,15.0000,17.8370];
n = [1,1,1,1,1,1,1,1];
r = [0,0,1,0,1,1,1,1];

nnum = sum(n);
%% 极大似然估计theta_e和对数似然函数
[theta_e, logL]= maximum_likelihood_estimates( x,n,r,nnum,datatype);

m_e = theta_e(1)
s_e = theta_e(2)*pi/sqrt(3.0)

%% Fisher 信息矩阵 其分量为i00，i01,i11
[i00,i01,i11] = fisher_information_matrix( x,n,r,nnum,theta_e,datatype)
%% Cramer-Rao下界，即

% Estimated Lower Bound of Variance of Mu，     varm
% Estimated Lower Bound of Variance of Sigma    vars
% Estimated CoVariance of Mu and Sigma          covms
[varm,vars,covms] = cramer_rao_bounds( x,n,r,nnum,theta_e,datatype);

varm = varm
vars = vars*pi*pi/3.0
covms = covms*pi/sqrt(3.0)
% 渐进方差置信域

%% 设置双侧置信度
cf_level = [ 0.5, 0.8, 0.9,0.95,0.98,0.99];
%% 概率为p发火点
p = [0.001,0.499];
%% 渐进方法置信域 m的置信区间 m_inter， s的置信域 s_inter Lp的置信域 lp_inter
[m_inter,s_inter,lp_inter] = asymptotic_confidence_intervals(x,n,r,nnum,datatype,cf_level,p);
m_inter = m_inter
s_inter = s_inter*pi/sqrt(3)
lp_inter = lp_inter

%% 画出渐进方法的联合置信域
plot_asymptotic_cr(theta_e,i00,i01,i11,m_inter,s_inter,datatype,cf_level)

%  方差函数方法 -- Langlie方法 置信区间求解方法验证。
[m_inter,s_inter,lp_inter,varm,vars,sigma] = variance_function_intervals(x,n,r,nnum,datatype,cf_level,p)
aaa = [m_inter,s_inter,lp_inter]
