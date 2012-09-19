function [m_inter,s_inter,lp_inter,varm,vars,sigma] = variance_function_intervals(x,n,r,nnum,datatype,cf_level,p)
%  ����fisher��Ϣ���� 
%   ���ݣ�x��n��r�� theta=��theta1��theta2��
%   ���룺
%             x :  �̼�ˮƽ
%             n :  ��ͬ�̼�ˮƽ�������
%             r :  ��ͬ�̼�ˮƽ��Ӧ����
%          nnum :  �̼�ˮƽ���� 
%         theta :  �����Ʋ���
%      datatype :  'norm'������̬�ֲ���'logistic'����logit�ֲ�
%      cf_level :  ���Ŷ�
%             p :  ��������
%   �����
%       m_inter :  ����m����������
%       s_inter :  ����s����������
%      lp_inter :  Lp����������
%          varm :  mu�ķ���
%          vars :  sigma�ķ���
%         sigma :  sigma����ƫ����
[theta_e,fval]= maximum_likelihood_estimates( x,n,r,nnum,datatype)
   
mu = theta_e(1);
sigma_e = theta_e(2);

%%%!!! ����������� Analysis of Sensitivity Tests ���µĹ�ʽ(6)������еļ��㷽����ͬ��
%beta = 1.6 /n^(2/3)     %% �����е�ȡ��
beta = 2/3               %% �����������Ϊ8ʱ��ȡ��

sigma = (1+beta)*sigma_e;

varm = 2.5*sigma^2 /nnum;
vars = 3.2*sigma^2 /nnum;
%vars = 1.2*sigma^2/n^(2/3)

cf = 1 - 0.5*(1 - cf_level);
cm1 = mu - sqrt(varm)*icdf(datatype,cf,0,1);
cm2 = mu + sqrt(varm)*icdf(datatype,cf,0,1);
cs1 = sigma - sqrt(vars)*icdf(datatype,cf,0,1);
cs2 = sigma + sqrt(vars)*icdf(datatype,cf,0,1);
m_inter = [cm1',cm2'];
s_inter = [cs1',cs2'];

lp_inter = [];
for i=1:length(p)
	zp = icdf(datatype,p(i),0,1);
	lp1 =  mu + zp*sigma - sqrt(varm + zp*zp *vars)*icdf(datatype,cf,0,1);
	lp2 =  mu + zp*sigma + sqrt(varm + zp*zp *vars)*icdf(datatype,cf,0,1);
	lp_inter = [lp_inter,lp1',lp2'];
end

end

