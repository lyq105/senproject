function ret =factorial( k )
%   k�Ľ׳�
  
if (k <= 1)
    ret = 1;
else
    ret =  factorial(k-1)*k;
end
end

