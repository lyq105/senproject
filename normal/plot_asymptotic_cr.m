function plot_asymptotic_cr(theta_0,i00,i01,i11,m_inter,s_inter,datatype,cf_level)
%% 画出渐进方法不同置信度的联合置信区域
   
nd = 200;
x0 = theta_0(1);
y0 = theta_0(2);
xmin = x0 - 1.2*(x0-min(min(m_inter)));
xmax = x0 + 1.2*(max(max(m_inter))-x0);

ymin = y0 - 1.2*(y0-min(min(s_inter)));
ymax = y0 + 1.2*(max(max(s_inter))-y0);

z = zeros(nd,nd);
x = linspace(xmin,xmax,nd);
y = linspace(ymin,ymax,nd);

for i = 1:nd
	for j = 1:nd
		z(i,j)=i00*(x(i) - x0)*(x(i) - x0) + 2*i01*(x(i)-x0)*(y(j)-y0) + i11*(y(j)-y0)*(y(j)-y0);
		z(i,j) = chi2cdf(z(i,j),2);
	end
end
d_cf = chi2inv(cf_level,1);
cf = chi2cdf(d_cf,2)

[C,h] = contour(x,y,z,cf);
clabel(C,cf);
hold on;
plot(x0,y0,'r^');
text_title = sprintf('Asymptotic Analysis Using Linear %s Response','Normal');
title(text_title);
xlabel('\mu');
ylabel('\sigma');

end

