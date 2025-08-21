%Script to perform source panel method on NACA 0012 form airfoil
clc; clear;
%========== Importing NACA 0012 Airfoil Profile and Generating Panels ==========%
U = 1; %free stream velocity
c = 1; %chord length
t_max = 0.12*c; %maximum thickness is 12% of chord length per naca 0012 airfoil profile
alphad = 11;
alpha = alphad*(pi/180);

[XB, YB, XC, YC, S, betaR, phiR] = loadFoil2(c, t_max, alphad);



%========== Geometric Integral Terms ==========%

[I J] = Calc_Iij(XC, YC, XB, YB, phiR, S); %I=normal integral term, J=tangent integral term

%========== Linear SoE Solution, Pressure Coefficient ==========%

[lambda, sum_lambda, V_s, Cp, NumPan, Gamma] = solvePanels(I, J, betaR, S, U);

disp(sum_lambda) %print the sum of lambda(j)*S(j) --- should be very close to zero
disp(Gamma)
%plot of discretized surface + control points

% indices = (V_s <= 0.0001); %finding the points of zero velocity, since we have set normal velocity = 0 already

figure; hold on; axis equal;
plot(XB,YB, 'b.', MarkerSize=7);
plot(XC, YC, 'r*');
plot(XB,YB,'k');
% plot(x_c(indices), y_c(indices), 'bo', MarkerSize=7, MarkerFaceColor='c')
title('Discretized Body Panels')
xlabel('X')
ylabel('Y')
legend('Panel Bounds', 'Control Points')



XB(end) = [];
half_x = floor(NumPan/2);
figure; hold on;
plot(XB(1:half_x), Cp(1:half_x), 'b');
plot(XB(half_x:end), Cp(half_x:end)); 
% plot(x_c, V_s, 'r--');
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphad), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Top Cp', 'Bottom Cp');






