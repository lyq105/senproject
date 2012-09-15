function plot_likelihood_ratio_cr(likelihood_func,theta_0,datatype,cf_level)
%% 画出渐进方法不同置信度的联合置信区域
   
nd = 20;
x0 = theta_0(1);
y0 = theta_0(2);
%xmin = x0 - 12*x0;
%xmax = x0 + 12*x0;
xmin = -50;
xmax = 20;

%ymin = y0 - 12*y0;
%ymax = y0 + 12*y0;

ymin = 0;
ymax = 200;

z = zeros(nd,nd);
x = linspace(xmin,xmax,nd);
y = linspace(ymin,ymax,nd);

for i = 1:nd
	for j = 1:nd
		z(i,j)= likelihood_func([x(i),sqrt(3.0)/pi*y(j)])/likelihood_func(theta_0);
		z(i,j) = chi2cdf(z(i,j),2)
	end
end
d_cf = chi2inv(cf_level,1);
cf = chi2cdf(d_cf,2)

%[C,h] = contour(x,y,z,cf);
[C,h] = contour(x,y,z,10);
%clabel(C,cf);
text_title = sprintf('Likelihood Ratio Analysis Using Linear %s Response',datatype);
title(text_title);
xlabel('\mu');
ylabel('\sigma');

end

