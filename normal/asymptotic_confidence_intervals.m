function [m_inter,s_inter,lp_inter] = asymptotic_confidence_intervals(x,n,r,nnum,datatype,cf_level,p)
%  ����fisher��Ϣ���� 
%   ���ݣ�x��n��r�� theta=��theta1��theta2��
%   ���룺
%             x :  �̼�ˮƽ
%             n :  ��ͬ�̼�ˮƽ�������
%             r :  ��ͬ�̼�ˮƽ��Ӧ����
%          nnum :  �̼�ˮƽ���� 
%         theta :  �����Ʋ���
%      datatype :  �ж��������� 0:��̬�ֲ���1��Logistic�ֲ�
%      cf_level :  ���Ŷ�
%             p :  ��������
%   �����
%       m_inter :  ����m����������
%       s_inter :  ����s����������
%      lp_inter :  Lp����������
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

