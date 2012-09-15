function plot_asymptotic_cr(theta_0,i00,i01,i11,m_inter,s_inter,datatype,cf_level)
%% 画出不同置信度的联合置信区域
   

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
		z(i,j)=i00*(x(i) - x0)*(x(i) - x0) + 2*i01*(x(i)-x0)*(y(j)*sqrt(3)/pi-y0) + i11*(y(j)*sqrt(3)/pi-y0)*(y(j)*sqrt(3)/pi-y0);
	end
end
d_cf = chi2inv(cf_level,1);
C = contour(x,y,z,d_cf);
clabel(C,chi2cdf(d_cf,2));
text_title = sprintf('Asymptotic Analysis Using Linear %s Response',datatype);
title(text_title);
xlabel('\mu');
ylabel('\sigma');
grid on;
%legend('10','10','10','10','10','10');

end

