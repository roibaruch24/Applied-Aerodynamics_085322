clear, close all
%% B.1.1
% Parameters
a = 0.5;
D = 2*a;
r0 = linspace(a, 2*a, 100);

% Calculate y0 for the Föppl lines
y0_pos = 0.5 * (r0 - (a^2 ./ r0));
y0_neg = -0.5 * (r0 - (a^2 ./ r0));

% Plot the cylinder
theta = linspace(0, 2*pi, 100);
x_cylinder = a * cos(theta);
y_cylinder = a * sin(theta);

% Vortex for Re = 20
x_0_r_20 = 0.83*D;
y_0_r_20 = 0.23*D;

% Vortex for Re = 1.43e5
x_0_r_1p4e5 = 0.81*D;
y_0_r_1p4e5 = 0.27*D;

figure;
hold on;
plot(x_cylinder, y_cylinder, 'k');
plot(r0, y0_pos, 'k', r0, y0_neg, 'k')
h1 = scatter(x_0_r_20,y_0_r_20,'filled','r');
scatter(x_0_r_20,-y_0_r_20,'filled','r');
h2 = scatter(x_0_r_1p4e5,y_0_r_1p4e5,'filled','b');
scatter(x_0_r_1p4e5,-y_0_r_1p4e5,'filled','b');
grid on
axis equal;
xlabel('x'); 
ylabel('y'); 
title('Cylinder with the theoretical Foppl line')

% Add legend for scatter points
legend([h1, h2], {'Re = 20', 'Re = 1.43\times 10^5'}, 'Location', 'best');
hold off;
%% B.1.2
x = linspace(-3*D,3*D);
y = linspace(-2*D,2*D);
[X,Y] = meshgrid(x,y);
z = X+1i*Y;
radius = abs(z);
inside_cylinder = radius < a;

w_Re_20 =z+0.25./z +1i*0.506*(log(z-0.83-1i*0.23)-log(z-0.83+1i*0.23)-log(z-0.25/(0.83-1i*0.23))+log(z-0.25/(0.83+1i*0.23)));
w_Re_1p4e5 = z+0.25./z +1i*0.495*(log(z-0.81-1i*0.27)-log(z-0.81+1i*0.27)-log(z-0.25/(0.81-1i*0.27))+log(z-0.25/(0.81+1i*0.27)));

w_Re_20(inside_cylinder) = NaN;
w_Re_1p4e5(inside_cylinder) = NaN;

figure;
subplot(1, 2, 1);
hold on;
plot(x_cylinder, y_cylinder, 'k');
contour(X, Y, imag(w_Re_20), 15);
title('Streamlines for Re = 20');
xlabel('$\frac{x}{D}$', 'Interpreter', 'latex');
ylabel('$\frac{y}{D}$', 'Interpreter', 'latex');
grid on;
axis equal;
hold off;


subplot(1, 2, 2);
hold on;
plot(x_cylinder, y_cylinder, 'k');
contour(X, Y, imag(w_Re_1p4e5), 15);
colormap('parula')
title('Streamlines for $Re = 1.4 \times 10^5$', 'Interpreter', 'latex');
xlabel('$\frac{x}{D}$', 'Interpreter', 'latex');
ylabel('$\frac{y}{D}$', 'Interpreter', 'latex');
grid on;
axis equal;
hold off;

figure;
hold on;
contour(X, Y, imag(w_Re_20), 15, 'r');
contour(X, Y, imag(w_Re_1p4e5), 15, 'b');
plot(x_cylinder, y_cylinder, 'k');
scatter(x_0_r_20,y_0_r_20,'filled','r');
scatter(x_0_r_20,-y_0_r_20,'filled','r');
scatter(x_0_r_1p4e5,y_0_r_1p4e5,'filled','b');
scatter(x_0_r_1p4e5,-y_0_r_1p4e5,'filled','b');
title('Streamlines for Different Reynolds Numbers');
xlabel('$\frac{x}{D}$', 'Interpreter', 'latex');
ylabel('$\frac{y}{D}$', 'Interpreter', 'latex');
legend('Re = 20', 'Re = 1.4 \times 10^5', 'Location', 'best');
grid on;
axis equal;
hold off;

%% B.2
load Cp_data.mat

% Parameters
theta = linspace(0,360,361);
z = a*cosd(theta)+1i*a*sind(theta);

gamma_Re_20 = 0.5061*D;
gamma_Re_1p4e5 = 0.495*D;

z0_Re_20 = 0.83*D + 0.23*D*i;  
z0_Re_1p4e5 = 0.81*D + 0.27*D*i;


% Calculate velocity around the cylinder
V_z_Re_20 = zeros(1, length(z));
V_z_Re_1p4e5 = zeros(1, length(z));

for i = 1:length(z)
    V_z_Re_20(i) = (1 - (a^2) / (z(i)^2)) + 1i * gamma_Re_20 * (1 / (z(i) - z0_Re_20) + 1 / (z(i) - a^2 / z0_Re_20) - 1 / (z(i) - conj(z0_Re_20)) - 1 / (z(i) - a^2 / conj(z0_Re_20)));
    V_z_Re_1p4e5(i) = (1 - (a^2) / (z(i)^2)) + 1i * gamma_Re_1p4e5 * (1 / (z(i) - z0_Re_1p4e5) + 1 / (z(i) - a^2 / z0_Re_1p4e5) - 1 / (z(i) - conj(z0_Re_1p4e5)) - 1 / (z(i) - a^2 / conj(z0_Re_1p4e5)));
end

Cp_foppl.Re_20 = 1 - abs(V_z_Re_20).^2;
Cp_foppl.Re_1p4e5 = 1 - abs(V_z_Re_1p4e5).^2;

% Plotting
figure;
hold on;
plot(data.theta{1},data.Cp{1});
plot(theta, Cp_foppl.Re_20);
xlabel('$\theta^\circ$ (deg)','Interpreter','latex');
ylabel('Cp');
title('Pressure Coefficient (Cp) around the Cylinder for Re = 36 - 107','Interpreter','latex');
legend('Literature Re = 36 - 107','Theortical Re = 36 - 107','Interpreter','latex');
grid on;
hold off;

figure;
hold on;
plot(data.theta{2},data.Cp{2});
plot(theta, Cp_foppl.Re_1p4e5);
xlabel('$\theta^\circ$ (deg)','Interpreter','latex');
ylabel('Cp');
title('Pressure Coefficient (Cp) around the Cylinder for $Re = 1.4 \times 10^5$','Interpreter','latex');
legend('Literature','Theortical','Interpreter','latex');
grid on;
hold off;



