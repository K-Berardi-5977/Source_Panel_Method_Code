function [lambda, sum_lambda, V_s, Cp, NumPan, Gamma, Vt_U] = solvePanels(I, J, beta, S, U)
NumPan = length(I(:,1)); %indexing vartiable
A = zeros(NumPan, NumPan); %variable to store Integral terms for linear system of equations

for i = 1:NumPan
    for j = 1:NumPan
        if (i == j)
            A(i,j) = pi; %normal velocity self-influence of the ith panel
        else
            A(i,j) = I(i,j);
        end
    end
end

%========== Free Stream Terms ==========%
b = zeros(NumPan,1); %variable to store free stream term for linear system of equations
for n = 1:NumPan
    b(n) = -U*2*pi*cos(beta(n));
end

%========== Final Source Strength Calculations ============%
lambda = A\b; %solving system of equations

%validation check, should be approximately zero
sum_lambda = zeros(NumPan,1);
for s = 1:length(lambda)
    sum_lambda(s) = lambda(s)*S(s);
end
sum_lambda=sum(sum_lambda);
%========== Tangent Velocity and Pressure Coefficient

V_s = zeros(NumPan,1);
Cp = zeros(NumPan, 1); 

for i = 1:NumPan
    source_terms = 0;
    for j = 1:NumPan
        source_terms = source_terms + (lambda(j)/(2*pi))*(J(i,j)); %accounts for sum of all source terms in the tangent velocity
    end
    Vt_U(i) = U*sin(beta(i));
    V_s(i) = Vt_U(i)+source_terms;
    Cp(i) = 1-(V_s(i)/U)^2;

    
end
Gamma = sum(V_s);
end

