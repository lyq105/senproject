function ret=likelihood_function( x,n,r,nnum,theta,datatype)
%   ���ָж����ݵ���Ȼ����
%   ���ݣ�x��n��r�� theta=��theta1��theta2��
%   x :  �̼�ˮƽ
%   n :  ��ͬ�̼�ˮƽ�������
%   r :  ��ͬ�̼�ˮƽ��Ӧ����
%   nnum : �̼�ˮƽ���� 
%   theta�������Ʋ���
%   datatype: �ж��������� 0:��̬�ֲ���1��Logistic�ֲ�

ret = 1;
for i =1:nnum
      z = (x(i) - theta(1)) / theta(2);
      ret = ret* combinatorialnumber(n(i),r(i)) * (cdf(datatype,z,0,1)^(r(i))) * ((1-cdf(datatype,z,0,1))^(n(i) - r(i)));  
end

end

