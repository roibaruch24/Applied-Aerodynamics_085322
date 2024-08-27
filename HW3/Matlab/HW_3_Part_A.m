function HW_3_Part_A(Input)

Input.L = 100;
Input.ni = 3.178e-5;
Input.x_0 = 1e-6;
Input.du_e_dx = 0;
Input.u_e = ((Input.Re_L) * Input.ni )/Input.L;
Input.Re_x_0 = (Input.u_e*Input.x_0)/Input.ni;
Input.theta_0 = (0.664*Input.x_0)/sqrt(Input.Re_x_0);
[Output_laminar] = laminar_boundary_layer(Input);
Input.H1_0 = 10.6;
[Output_turbulent] = turbulent_boundary_layer(Input);

% Plotting 
figure('Name', 'theta vs x')
plot(Output_laminar.x_l/Input.L, Output_laminar.theta_Thwaites);
xlabel('$\frac{x}{L}$', 'Interpreter', 'latex');
ylabel('$\theta$', 'Interpreter', 'latex');
grid on
hold on
plot(Output_laminar.x_l/Input.L, Output_laminar.theta_Blasius)
title(['$\theta$ vs x for Re = ' num2str(Input.Re_L)], 'Interpreter', 'latex')
legend('Thwaites', 'Blasius')

figure('Name', 'delta star vs x')
plot(Output_laminar.x_l/Input.L, Output_laminar.delta_star_Thwaites);
xlabel('$\frac{x}{L}$', 'Interpreter', 'latex');
ylabel('$\delta^*$', 'Interpreter', 'latex');
grid on
hold on
plot(Output_laminar.x_l/Input.L, Output_laminar.delta_star_Blasius)
title(['$\delta^*$ vs x for Re = ' num2str(Input.Re_L)], 'Interpreter', 'latex')
legend('Thwaites', 'Blasius')

figure('Name', 'cf vs x')
semilogy(Output_laminar.x_l/Input.L, Output_laminar.cf_Thwaites);
xlabel('$\frac{x}{L}$', 'Interpreter', 'latex');
ylabel('$c_f$', 'Interpreter', 'latex');
grid on
hold on
semilogy(Output_laminar.x_l/Input.L, Output_laminar.cf_Blasius)
title(['$c_f$ vs x for Re = ' num2str(Input.Re_L)], 'Interpreter', 'latex')
legend('Thwaites', 'Blasius')

% Plotting
figure('Name', 'H vs x')
plot(Output_laminar.x_l/Input.L, Output_laminar.H_Thwaites);
xlabel('$\frac{x}{L}$', 'Interpreter', 'latex');
ylabel('H', 'Interpreter', 'latex');
grid on
hold on
plot(Output_laminar.x_l/Input.L, Output_laminar.H_Blasius)
title(['H vs x for Re = ' num2str(Input.Re_L)], 'Interpreter', 'latex')
legend('Thwaites', 'Blasius')



figure('Name','theta turbulent vs x')
plot(Output_turbulent.x_t/Input.L, Output_turbulent.theta_head);
xlabel('$\frac{x}{L}$');
ylabel('$\theta$');
title(['$\theta$ vs x for Re = ' num2str(Input.Re_L)])
grid on
hold on
plot(Output_turbulent.x_t/Input.L, Output_turbulent.theta_Schlichting);
legend('Head', 'Schlichting')

figure('Name','H turbulent vs x')
semilogx(Output_turbulent.x_t/Input.L, Output_turbulent.H_head);
xlabel('$\frac{x}{L}$');
ylabel('H');
title(['H vs x for Re = ' num2str(Input.Re_L)])
grid on
hold on
semilogx(Output_turbulent.x_t/Input.L, Output_turbulent.H_Schlichting);
legend('Head', 'Schlichting')

figure('Name','cf turbulent vs x')
semilogx(Output_turbulent.x_t/Input.L, Output_turbulent.cf_head);
xlabel('$\frac{x}{L}$');
ylabel('cf');
title(['$c_f$ vs x for Re = ' num2str(Input.Re_L)])
grid on
hold on
semilogx(Output_turbulent.x_t/Input.L, Output_turbulent.cf_Schlichting)
legend('Head', 'Schlichting')

figure('Name','delta star turbulent vs x')
plot(Output_turbulent.x_t/Input.L, Output_turbulent.delta_star_head);
xlabel('$\frac{x}{L}$');
ylabel('$\delta^*$');
title(['$\delta^*$ vs x for Re = ' num2str(Input.Re_L)])
grid on
hold on
plot(Output_turbulent.x_t/Input.L, Output_turbulent.delta_star_Schlichting);
legend('Head', 'Schlichting')
end