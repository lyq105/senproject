function [theta_e,fval]=maximum_likelihood_estimates( x,n,r,nnum,datatype)
% ���㼫����Ȼ����
%   ���ݣ�x��n��r�� theta=��theta1��theta2��
%   x :  �̼�ˮƽ
%   n :  ��ͬ�̼�ˮƽ�������
%   r :  ��ͬ�̼�ˮƽ��Ӧ����
%   nnum : �̼�ˮƽ���� 
%   theta�������Ʋ���
%   datatype: �ж��������� 'norm':��̬�ֲ���'logistic'��Logistic�ֲ�

[theta_e,fval] = fminsearch(@(theta)-log(likelihood_function(x,n,r,nnum,theta,datatype)),[mean(x),var(x)]);
fval = - fval
end 
