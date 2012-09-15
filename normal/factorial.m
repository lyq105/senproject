function ret =factorial( k )
%   kµÄ½×³Ë
  
if (k <= 1)
    ret = 1;
else
    ret =  factorial(k-1)*k;
end
end

