function [Input] = Panel_code_init(Input)
% This part defines the NACA airfoil and the Input parametrs, only works for the 4 digit series
Input.N = 150;
Input.Chord = 1; % m
Input.V_inf  = 1; % m/s


N = (Input.N/2)+1;
X_airfoil = zeros(1,N);
tmax = str2double(Input.NACA(3:4))/100; %the maximum airfoil x thickness 
for i = 1:(N)
    X_airfoil(i) = 0.5*(1-cos((i-1)*pi/(N-1)));
end
switch Input.TE_type
    case 'blunt'
        Y_airfoil = 5 * tmax * (0.2969 * sqrt(X_airfoil) - ...
            0.1260 *  X_airfoil - ...
            0.3516 *  X_airfoil.^2 + ...
            0.2843 *  X_airfoil.^3 - ...
            0.1015 *  X_airfoil.^4);
    case 'sharp'
        Y_airfoil = 5 * tmax * (0.2969 * sqrt(X_airfoil) - ...
            0.1260 *  X_airfoil - ...
            0.3516 *  X_airfoil.^2 + ...
            0.2843 *  X_airfoil.^3 - ...
            0.1036 *  X_airfoil.^4);
end

% Mean Camber
X_camber = zeros(1,N);
Y_camber = zeros(1,N);
m = str2double(Input.NACA(1))/100; % maximum camber
p = str2double(Input.NACA(2))/10; % maximum_camber_location
for i = 1:(N)
    X_camber(i) = 0.5*(1-cos((i-1)*pi/(N-1)));
end
for i = 1:length(X_camber)
    if (X_camber(i) <= p)
        Y_camber(i) = (m*X_camber(i)*(2*p-X_camber(i)))/p^2;
    else
       Y_camber(i) = (m*(1-2*p+2*p*X_camber(i)-X_camber(i)^2))/((1-p)^2);
    end 
end
%NACA upper and lower vectors + ploting the airfoil
theta = atan(gradient(Y_camber));
X_u = X_airfoil - Y_airfoil.*sin(theta);
X_l = X_airfoil + Y_airfoil.*sin(theta);
Y_u = Y_camber + Y_airfoil.*cos(theta);
Y_l = Y_camber - Y_airfoil.*cos(theta);

% Input X, Y vectors  
Input.X_vec = [fliplr(X_l) X_u(2:end)];
Input.Y_vec = [fliplr(Y_l) Y_u(2:end)];



end