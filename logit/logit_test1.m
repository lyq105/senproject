%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     logistic�ֲ����ݼ��㷽����֤
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ����
datatype = 'logistic';

x = [3.00000,6.07538,7.41994,9.00000,10.9945,12.0000,15.0000,17.8370];
n = [1,1,1,1,1,1,1,1];
r = [0,0,1,0,1,1,1,1];

nnum = sum(n);

%% ������Ȼ����
[theta_e, logL]= maximum_likelihood_estimates( x,n,r,nnum,datatype)
%% Fisher ��Ϣ����
[i00,i01,i11] = fisher_information_matrix( x,n,r,nnum,theta_e,datatype)
%% Cramer-Rao�½� 
[varm,vars,covms] = cramer_rao_bounds( x,n,r,nnum,theta_e,datatype)

%% ��������������

%% ���Ŷ�
cf_level = [ 0.5, 0.8, 0.9,0.95,0.98,0.99];
%% ����p
p = [0.001,0.999];
%% ��������������
[m_inter,s_inter,lp_inter] = asymptotic_confidence_intervals(x,n,r,nnum,datatype,cf_level,p)
%% ��������������
plot_asymptotic_cr(theta_e,i00,i01,i11,m_inter,s_inter,datatype,cf_level)

