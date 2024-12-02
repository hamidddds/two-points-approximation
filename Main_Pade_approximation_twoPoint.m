% Clear command window
clc;
clear all;


syms u;  % Define u as a symbolic variable
% Calculate the summation
main_func = 0;
for n = 0:7
    term = ((-1)^n / factorial(n)) * (2)^(n/2) * gamma(1 + 0.5 * n) * u^n;
    main_func = main_func + term;
end
% Display the summation result
main_func=vpa(main_func,2);
[expanded_coeffs, expanded_terms] = coeffs(main_func, u);
ZeroSide_coeff=flip(expanded_coeffs);

Inf_Side_Coeff=[0,0,1];


% Example usage
M_denominator_count = 5;    % Number of coefficients for the denominator

N_Zero_side_coeffs = length(ZeroSide_coeff)-1               % Maximum power of u in the expanded series zero side
N_Inf_sideCoeff = length(Inf_Side_Coeff)-1     % Maximum power of u in the expanded series for comparison in inifinity


% Parameters for the summation
sigma = 1;




c = sym('c', [1 N_Zero_side_coeffs+1]); % Create a symbolic array with n+1 symbols

b = sym('b', [1 N_Zero_side_coeffs+1]); % Create a symbolic array with n+1 symbols
a = sym('a', [1 N_Zero_side_coeffs+1]); % Create a symbolic array with n+1 symbols
d = sym('d', [1 N_Inf_sideCoeff+1]); % Create a symbolic array with n+1 symbols
eq = sym([]);
eq1 = sym([]);





for j = 0:N_Zero_side_coeffs
    temp_a = 0;
    for i=0:j
        temp_a = temp_a + b(i+1)*c(j-i+1); 
    end
    eq(j+1)=temp_a-a(j+1)==0;
end


for j = 0:N_Zero_side_coeffs
    for i = M_denominator_count+1:N_Zero_side_coeffs 
        eq(j+1)=subs(eq(j+1), [b(i+1),b(1)], [0,1]); % Substitute b(i) with 0 in a(j)
    end
end

for j = 0:N_Zero_side_coeffs
    for i = M_denominator_count:N_Zero_side_coeffs 
        eq(j+1)=subs(eq(j+1), a(i+1), 0); % Substitute b(i) with 0 in a(j)
    end
end



k=1;
for l = 0:N_Inf_sideCoeff-1
    temp_a = 0;
    for i=0:l
        temp_a=temp_a+d(i+1+1)*b(M_denominator_count-l+1+i); 
    end
    eq1(k)=temp_a-a(M_denominator_count-l-1+1)==0;    
    k=k+1;
end



k=1;
for i=length(eq)+1:length(eq)+length(eq1)
    eq(i)=eq1(k);
    k=k+1;
end





for i= 1:length(eq)
    eq(i) = subs(eq(i), c, ZeroSide_coeff);
    eq(i) = subs(eq(i), d, Inf_Side_Coeff);
end

symbols_a = symvar(eq);
solutions = solve(eq, symbols_a)


% Generate the numerator and denominator polynomials
[numerator, denominator] = create_polynomial(M_denominator_count, M_denominator_count,u);

% Display the rational function
rational_func = numerator / denominator





% Initialize an empty array to store the values
values = [];

% Loop through all fields in the structure
fields = fieldnames(solutions);
for i = 1:length(fields)
    % Extract each field value and append it to the values array
    values = [values, solutions.(fields{i})];
end


values = vpa(values,5)

symbols_a
rational_func

result = subs(rational_func, symbols_a, values);

xxx=subs(result, u, 5)
x = 5;
xxxx = (0.01*x^3 + 0.11*x^2 + 0.47*x + 1) / (0.01*x^5 + 0.11*x^4 + 0.51*x^3 + 1.27*x^2 + 1.72*x + 1)



function [numerator, denominator] = create_polynomial(numerator_count, denominator_count,u)
    % Define symbolic variable      

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


