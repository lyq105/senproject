function ret = combinatorialnumber( a,b )
%   �����C_a^b

ret = (factorial(a) / factorial(b) / factorial(a-b));
end

