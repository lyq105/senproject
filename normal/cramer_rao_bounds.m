function [varm,vars,covms] = cramer_rao_bounds( x,n,r,nnum,theta,datatype)
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

[i00,i01,i11] = fisher_information_matrix(x,n,r,nnum,theta,datatype);

deti  =  i00*i11- i01*i01;
varm  =  i11/deti;	
vars  =  i00/deti;	
covms = -i01/deti;

end

