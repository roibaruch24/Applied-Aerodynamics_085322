function [X_camber, Y_camber] = Mean_Camber(Input)
%
N = (Input.N/2)+1;
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
end