function [X_u, X_l, Y_u, Y_l] = NACA_vectors(X_airfoil, Y_airfoil, X_camber, Y_camber, Input)
%
theta = atan(gradient(Y_camber));
X_u = X_airfoil - Y_airfoil.*sin(theta);
X_l = X_airfoil + Y_airfoil.*sin(theta);
Y_u = Y_camber + Y_airfoil.*cos(theta);
Y_l = Y_camber - Y_airfoil.*cos(theta);
figure;
hold on;
grid on;
axis equal;
scatter(X_l, Y_l)
scatter(X_u, Y_u)
title(['NACA ' Input.NACA])
xlabel('$\frac{x}{D}$')
ylabel('$\frac{y}{D}$')
plot(X_camber, Y_camber)
hold off;
legend('Upper','Lower','Mean')
end