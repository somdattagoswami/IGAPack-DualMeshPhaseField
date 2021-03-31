function [y] = exact_sol_phi(x,Fract)

y = zeros(size(x));
y = exp(-abs(x-0.0)/Fract.constl);

end