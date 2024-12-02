% function [Estimated_func] = TwoPointPade_main(u,DegreeDenominator, ZeroCoefficient, InfCoefficient,Coefficient_J,Coefficient_K)
syms u;
DegreeDenominator = 2
ZeroCoefficient = 3
InfCoefficient = [1,-3,6.5]
InfCoefficient = [0]
Coefficient_J=[1,-3,6.5]
Coefficient_K=[0]

% the form of approximation is  f(x) = (a*x^(M-1) + b*x^(M-2)+ ... + c) /
% (d*x^M + e*x^M-1 + ...+ 1) . with this function we approximat an arbitary
% function with 2 points. one point must be define in infinity and the
% second point is define near zero.

J_zero = ZeroCoefficient;               % Maximum power of u in the expanded series for comparison
K_infinity = InfCoefficient;     % Maximum power of u in the expanded series for comparison in inifinity

c = sym('c', [1 ZeroCoefficient+1]); % Create a symbolic array with n+1 symbols
b = sym('b', [1 ZeroCoefficient+1]); % Create a symbolic array with n+1 symbols
a = sym('a', [1 ZeroCoefficient+1]); % Create a symbolic array with n+1 symbols
d = sym('d', [1 K_infinity+1]); % Create a symbolic array with n+1 symbols
eq = sym([]);
eq1 = sym([]);


for j = 0:J_zero
    temp_a = 0;
    for i=0:j
        temp_a = temp_a + b(i+1)*c(j-i+1); 
    end
    eq(j+1)=temp_a-a(j+1)==0;
end


for j = 0:J_zero
    for i = DegreeDenominator+1:J_zero 
        eq(j+1)=subs(eq(j+1), [b(i+1),b(1)], [0,1]); % Substitute b(i) with 0 in a(j)
    end
end

for j = 0:J_zero
    for i = DegreeDenominator:J_zero 
        eq(j+1)=subs(eq(j+1), a(i+1), 0); % Substitute b(i) with 0 in a(j)
    end
end



k=1;
for l = 0:K_infinity-1
    temp_a = 0;
    for i=0:l
        temp_a=temp_a+d(i+1+1)*b(DegreeDenominator-l+1+i); 
    end
    eq1(k)=temp_a-a(DegreeDenominator-l-1+1)==0;    
    k=k+1;
end



k=1;
for i=length(eq)+1:length(eq)+length(eq1)
    eq(i)=eq1(k);
    k=k+1;
end

for i= 1:length(eq)
    eq(i) = subs(eq(i), c, Coefficient_J);
    eq(i) = subs(eq(i), d, Coefficient_K);
end

symbols_a = symvar(eq);
solutions = solve(eq, symbols_a);


% Generate the numerator and denominator polynomials

[numerator, denominator] = create_polynomial(DegreeDenominator, DegreeDenominator,u);

% Display the rational function
rational_func = numerator / denominator;


% Initialize an empty array to store the values
values = [];

% Loop through all fields in the structure
fields = fieldnames(solutions);
for i = 1:length(fields)
    % Extract each field value and append it to the values array
    values = [values, solutions.(fields{i})];
end



Estimated_func = subs(rational_func, symbols_a, values);

xxx=subs(Estimated_func, u, 5)
x = 5;
xxxx = (0.01*x^3 + 0.11*x^2 + 0.47*x + 1) / (0.01*x^5 + 0.11*x^4 + 0.51*x^3 + 1.27*x^2 + 1.72*x + 1)



function [numerator, denominator] = create_polynomial(numerator_count, denominator_count,u)
    
    % Create symbolic coefficients for the numerator starting from a_0
    numerator_coeffs = sym('a', [1, numerator_count]);
   
    % Create the numerator polynomial
    numerator = sum(numerator_coeffs .* (u.^(0:numerator_count-1)));

    % Define the coefficients for the denominator starting from b_1
    denominator_coeffs = sym('b', [1, denominator_count]);
    for i = 1:denominator_count
        denominator_coeffs(i) = sym(sprintf('b%d', i + 1));
    end

    % Create the denominator polynomial (assuming b0 is 1)
    denominator = 1 + sum(denominator_coeffs .* u.^(1:denominator_count));

end



% end

