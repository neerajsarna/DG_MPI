function f = HermiteHalfSpace(m,n)

if m == n
    f = 0.5;
end

if mod((m-n),2)==0 && m~=n
    f = 0;
end

if mod((m-n),2) ~= 0
  fact1 = sqrt(factorial(m) / (factorial(n) * 2 * pi));  
  if m > n
      fact2 = factd(m-n-2) * (-1)^((m-n-1)/2) / factorial(m-n);
  else
      fact2 = 0;
  end
  
  
  fact3 = compute_sum_gamma(m,n);
  
  f = (fact2 + fact3) * fact1; 
  
end

end

function f = compute_sum_gamma(m,n)
f = 0 ;

for k = 1 : n
    f = f + pi * sqrt(2)^(1-2 * k + m + n)/(gamma((k-m)/2)*gamma(2-k+m)*gamma((1+k-n)/2));
end
end