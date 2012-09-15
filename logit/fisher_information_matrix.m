function [i00,i01,i11] = fisher_information_matrix( x,n,r,nnum,theta,datatype)
%  ����fisher��Ϣ���� 
%   ���ݣ�x��n��r�� theta=��theta1��theta2��
%   ���룺
%       x :  �̼�ˮƽ
%       n :  ��ͬ�̼�ˮƽ�������
%       r :  ��ͬ�̼�ˮƽ��Ӧ����
%       nnum : �̼�ˮƽ���� 
%       theta�������Ʋ���
%       datatype: �ж��������� 0:��̬�ֲ���1��Logistic�ֲ�
%   �����
%       ��Ϣ����I = [i00,i01;i01,i11]
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

