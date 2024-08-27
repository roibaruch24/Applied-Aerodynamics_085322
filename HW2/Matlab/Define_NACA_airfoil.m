function [X_airfoil, Y_airfoil] = Define_NACA_airfoil(Input)
% Defines a given Input.NACA 4-digit 
N = (Input.N/2)+1;
X_airfoil = zeros(1,N);
tmax = str2double(Input.NACA(3:4))/100; %the maximum airfoil x thickness 
for i = 1:(N)
    X_airfoil(i) = 0.5*(1-cos((i-1)*pi/(N-1)));
end
Y_airfoil = 5 * tmax * (0.2969 * sqrt(X_airfoil) - ...
            0.1260 *  X_airfoil - ...
            0.3516 *  X_airfoil.^2 + ...
            0.2843 *  X_airfoil.^3 - ...
            0.1015 *  X_airfoil.^4); 
end