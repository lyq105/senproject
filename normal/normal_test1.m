%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     正态分布数据计算方法验证
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 数据
datatype = 'norm';
x = [10.0000,11.2500,11.8750,12.1875,12.5000,12.8812,15.0000,20.0000];
n = [ 1,1,1,1,1,1,1,1];
r = [0,0,0,1,1,0,1,1];
nnum = sum(n);

%% 极大似然估计
[theta_e, logL]= maximum_likelihood_estimates( x,n,r,nnum,datatype)
%% Fisher 信息矩阵
[i00,i01,i11] = fisher_information_matrix( x,n,r,nnum,theta_e,datatype)
%% Cramer-Rao下界 
[varm,vars,covms] = cramer_rao_bounds( x,n,r,nnum,theta_e,datatype)

%% 渐进方差置信域

%% 置信度
cf_level = [ 0.5, 0.8, 0.9,0.95,0.98,0.99];
%% 概率p
p = [0.001,0.999];
%% 渐进方法置信域
[m_inter,s_inter,lp_inter] = asymptotic_confidence_intervals(x,n,r,nnum,datatype,cf_level,p)
%% 画出联合置信域
plot_asymptotic_cr(x,n,r,nnum,datatype,cf_level,p)


