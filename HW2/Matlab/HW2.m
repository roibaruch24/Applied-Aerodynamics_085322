clear, close all
%% Main
alpha_vec = 0:1:10;
N_vec = [10:20:100 200 400];

% Preallocate arrays
Cp_Hess_Smith = zeros(length(N_vec), length(alpha_vec), max(N_vec));
CL_Hess_Smith = zeros(length(N_vec), length(alpha_vec));
CD_Hess_Smith = zeros(length(N_vec), length(alpha_vec));
CM_Hess_Smith = zeros(length(N_vec), length(alpha_vec));
CL_kutta_Hess_Smith = zeros(length(N_vec), length(alpha_vec));
CM_kutta_Hess_Smith = zeros(length(N_vec), length(alpha_vec));

Cp_vpm = zeros(length(N_vec), length(alpha_vec), max(N_vec)); 
CL_vpm = zeros(length(N_vec), length(alpha_vec));
CD_vpm = zeros(length(N_vec), length(alpha_vec));
CM_vpm = zeros(length(N_vec), length(alpha_vec));
CL_kutta_vpm = zeros(length(N_vec), length(alpha_vec));
CM_kutta_vpm = zeros(length(N_vec), length(alpha_vec));

for i = 1:length(N_vec)
    for j = 1:length(alpha_vec)
        Input.N = N_vec(i);
        Input.TE_type = 'blunt'; % can be either blunt or sharp
        [Input] = Panel_code_init(Input);
        Input.alpha = alpha_vec(j);
        [Cp_Hess_Smith(i, j, 1:Input.N), CL_Hess_Smith(i, j), CD_Hess_Smith(i, j), CM_Hess_Smith(i, j), CL_kutta_Hess_Smith(i, j),CM_kutta_Hess_Smith(i,j), X_mid_Hess_smith, U_tanget(i, j)] = Hess_Smith(Input);
        [Cp_vpm(i, j, 1:Input.N), CL_vpm(i, j), CD_vpm(i, j), CM_vpm(i, j),CL_kutta_vpm(i,j), CM_kutta_vpm(i,j), X_mid_vpm] = Vortex_panel_method(Input);
    end
end

%% Plotting CL, CD, CM vs alpha

figure('Name', 'C_L vortex panel method');
hold on
for i = 1:length(N_vec)
    plot(alpha_vec, CL_vpm(i,:),'DisplayName', ['N = ' num2str(N_vec(i))])
end
grid on
legend show
title('$C_L$ vs angle of attack', 'Interpreter', 'latex')
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$C_L$', 'Interpreter', 'latex')

figure('Name', 'C_D vortex panel method');
hold on
for i = 1:length(N_vec)
    plot(alpha_vec, CD_vpm(i,:),'DisplayName', ['N = ' num2str(N_vec(i))])
end
grid on
legend show
title('$C_D$ vs angle of attack', 'Interpreter', 'latex')
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$C_D$', 'Interpreter', 'latex')

figure('Name', 'C_M vortex panel method');
hold on
for i = 1:length(N_vec)
    plot(alpha_vec, CM_vpm(i,:),'DisplayName', ['N = ' num2str(N_vec(i))])
end
grid on
legend show
title('$C_M$ vs angle of attack', 'Interpreter', 'latex')
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$C_M$', 'Interpreter', 'latex')

%% Plot CL, CD, CM vs N

figure('Name', 'Coefficients Convergence Hess-Smith');
hold on
plot(N_vec, CL_Hess_Smith(:,5))
plot(N_vec, CD_Hess_Smith(:,5))
plot(N_vec, CM_Hess_Smith(:,5))
plot(N_vec, CL_kutta_Hess_Smith(:,5))
plot(N_vec, CM_kutta_Hess_Smith(:,5))
title('$C_L$, $C_D$, $C_M$ vs number of panels Hess Smith', 'Interpreter', 'latex')
xlabel('N', 'Interpreter', 'latex')
grid on
legend ('$C_L$', '$C_D$', '$C_M$','$C_L$ kutta','$C_M$ kutta', 'Interpreter', 'latex')

figure('Name', 'Coefficients Convergence Vortex Panel Method');
hold on
plot(N_vec, CL_vpm(:,5))
plot(N_vec, CD_vpm(:,5))
plot(N_vec, CM_vpm(:,5))
plot(N_vec, CL_kutta_vpm(:,5))
plot(N_vec, CM_kutta_vpm(:,5))
title('$C_L$, $C_D$, $C_M$ vs number of panels Vortex panel method', 'Interpreter', 'latex')
xlabel('N', 'Interpreter', 'latex')
grid on
legend ('$C_L$', '$C_D$', '$C_M$','$C_L$ kutta','$C_M$ kutta', 'Interpreter', 'latex')

%% Comparison of Methods
figure('Name', 'Coefficients Comparison');
subplot(3,1,1)
hold on
plot(N_vec, CL_vpm(:,5))
plot(N_vec, CL_Hess_Smith(:,5))
grid on
ylabel('$C_L$', 'Interpreter', 'latex')
xlabel('N', 'Interpreter', 'latex')
title('$C_L$ Convergence for Different Panel Methods', 'Interpreter', 'latex')
legend('Vortex panel method','Hess Smith', 'Interpreter', 'latex')

subplot(3,1,2)
hold on
plot(N_vec, CD_vpm(:,5))
plot(N_vec, CD_Hess_Smith(:,5))
grid on
ylabel('$C_D$', 'Interpreter', 'latex')
xlabel('N', 'Interpreter', 'latex')
title('$C_D$ Convergence for Different Panel Methods', 'Interpreter', 'latex')
legend('Vortex panel method','Hess Smith', 'Interpreter', 'latex')

subplot(3,1,3)
hold on
plot(N_vec, CM_vpm(:,5))
plot(N_vec, CM_Hess_Smith(:,5))
grid on
ylabel('$C_M$', 'Interpreter', 'latex')
xlabel('N', 'Interpreter', 'latex')
title('$C_M$ Convergence for Different Panel Methods', 'Interpreter', 'latex')
legend('Vortex panel method','Hess Smith', 'Interpreter', 'latex')

%%
clear
Input.TE_type = 'blunt'; % can be either blunt or sharp
Input.N = 150;
[Input] = Panel_code_init(Input);
Input.alpha = 10;
[Cp_Hess_Smith, CL_Hess_Smith, CD_Hess_Smith, CM_Hess_Smith, CL_kutta_Hess_Smith,CM_kutta_Hess_Smith, X_mid_Hess_smith] = Hess_Smith(Input);
[Cp_vpm, CL_vpm, CD_vpm, CM_vpm,CL_kutta_vpm, CM_kutta_vpm, X_mid_vpm] = Vortex_panel_method(Input);

figure('Name', 'C_P for Both Methods');
hold on
plot(X_mid_Hess_smith, -Cp_Hess_Smith)
plot(X_mid_vpm, -Cp_vpm)
grid on
xlabel('$\frac{x}{c}$', 'Interpreter', 'latex')
ylabel('$-C_P$', 'Interpreter', 'latex')
title('$C_P$ at $\alpha = 10$', 'Interpreter', 'latex')
legend('Hess Smith','Vortex panel method', 'Interpreter', 'latex')

%%
Input.N = 150;
Input.TE_type = 'blunt';
[Input] = Panel_code_init(Input);
Airfoil_Plotter(Input)

Input.N = 150;
Input.TE_type = 'sharp';
[Input] = Panel_code_init(Input);
Airfoil_Plotter(Input)

%% Comparison of TE types
clear
figure('Name', 'TE types Comparison');
Input.N = 300;
alpha_vec = 0:1:10;
Input.TE_type = 'blunt';
[Input] = Panel_code_init(Input);
for i = 1:length(alpha_vec)
    Input.alpha = alpha_vec(i);
    Input.TE_type = 'blunt';
    [Input] = Panel_code_init(Input);
    [Cp_Hess_Smith_blunt(i,:), CL_Hess_Smith_blunt(i), CD_Hess_Smith_blunt(i), CM_Hess_Smith_blunt(i), CL_kutta_Hess_Smith_blunt(i), CM_kutta_Hess_Smith_blunt(i), X_mid_Hess_smith_blunt] = Hess_Smith(Input);
    Input.TE_type = 'sharp';
    [Input] = Panel_code_init(Input);
    [Cp_Hess_Smith_sharp(i,:), CL_Hess_Smith_sharp(i), CD_Hess_Smith_sharp(i), CM_Hess_Smith_sharp(i), CL_kutta_Hess_Smith_sharp(i),CM_kutta_Hess_Smith_sharp(i), X_mid_Hess_smith_sharp] = Hess_Smith(Input);
end

subplot(3,1,1)
plot(alpha_vec, CL_Hess_Smith_blunt)
hold on
plot(alpha_vec, CL_Hess_Smith_sharp)
grid on
ylabel('$C_L$', 'Interpreter', 'latex')
xlabel('$\alpha$', 'Interpreter', 'latex')
title('$C_L$ Convergence for Different Panel Methods', 'Interpreter', 'latex')
legend('Blunt','Sharp', 'Interpreter', 'latex')

subplot(3,1,2)
plot(alpha_vec, CD_Hess_Smith_blunt)
hold on
plot(alpha_vec, CD_Hess_Smith_sharp)
grid on
ylabel('$C_D$', 'Interpreter', 'latex')
xlabel('$\alpha$', 'Interpreter', 'latex')
title('$C_D$ Convergence for Different Panel Methods', 'Interpreter', 'latex')
legend('Blunt','Sharp', 'Interpreter', 'latex')

subplot(3,1,3)
plot(alpha_vec, CM_Hess_Smith_blunt)
hold on
plot(alpha_vec, CM_Hess_Smith_sharp)
grid on
ylabel('$C_M$', 'Interpreter', 'latex')
xlabel('$\alpha$', 'Interpreter', 'latex')
title('$C_M$ Convergence for Different Panel Methods', 'Interpreter', 'latex')
legend('Blunt','Sharp', 'Interpreter', 'latex')

%%
clear
figure('Name','Thin airfoil theory')
alpha_vec = 0:1:15;
c_l = 2*pi.*alpha_vec + 0.227;
c_m = zeros(1,length(alpha_vec)) -0.053;
plot(alpha_vec, c_l)
hold on
plot(alpha_vec, c_m)
ylim([-10 100])
grid on
xlabel('$\alpha$', 'Interpreter', 'latex')
title('$C_l$ and $C_m$ for thin airfoil theory', 'Interpreter', 'latex')
legend('$C_l$','$C_m$', 'Interpreter', 'latex')

%% 
clear
alpha = 10*pi/180; 
A_n_vec(1) = alpha -0.0046;
i=1;
while abs(A_n_vec(end)) >= 1*10^-3
    syms x
    exp1 = cos(i*x)*(0.125*cos(x)-0.025);
    exp2 = cos(i*x)*(0.0555*cos(x)-0.011);
    I1 = int(exp1,x,[0 acos(1/5)]);
    I2 = int(exp2,x,[acos(1/5) pi]);
    I1 = double(I1);
    I2 = double(I2);
    A_n_vec(end+1) = (2/pi)*(I1+I2);
    i=i+1;
end

phi = (10:1:170)/57.3;
gamma = zeros(1,length(phi)+1);

for i = 1:length(phi)
    a =phi(i);
    b= rad2deg(phi(i));
    gamma(i) = 2*A_n_vec(1)*(1+cos(phi(i)))/sin(phi(i));
    for j = 1:length(A_n_vec)-1
        gamma(i) = (gamma(i)+2*A_n_vec(j+1)*sin(j*phi(i)));
    end
end

Cp_thin_airfoil = 2*gamma;
phi = [phi,pi];
x = 0.5*(1-cos(phi));


figure('Name','Cp thin airfoil')
plot(x,Cp_thin_airfoil)
title('$C_p$ using thin airfoil theory','Interpreter','latex')
xlabel('$\frac{x}{c}$','Interpreter','latex')
ylabel('$\Delta C_p$','Interpreter','latex')
grid on


Input.TE_type = 'blunt'; % can be either blunt or sharp
Input.N = 150;
Input.alpha = 10;
[Input] = Panel_code_init(Input);

figure('Name','Cp thin airfoil and panels')
plot(x,Cp_thin_airfoil)
hold on
[Cp_Hess_Smith_blunt, ~, ~, ~, ~, ~, X_mid_Hess_smith_blunt] = Hess_Smith(Input);
Input.TE_type = 'sharp';
[Input] = Panel_code_init(Input);
[Cp_Hess_Smith_sharp, ~, ~, ~, ~, ~, X_mid_Hess_smith_sharp] = Hess_Smith(Input);
plot(X_mid_Hess_smith_blunt, -Cp_Hess_Smith_blunt)
plot(X_mid_Hess_smith_sharp, -Cp_Hess_Smith_sharp)
title('$C_p$ using thin airfoil theory','Interpreter','latex')
xlabel('$\frac{x}{c}$','Interpreter','latex')
ylabel('$\Delta C_p$','Interpreter','latex')
grid on
hold off
legend('Thin airfoil','Blunt TE','Sharp TE')


