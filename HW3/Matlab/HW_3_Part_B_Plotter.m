function HW_3_Part_B_Plotter(Output, Reynolds_list, alpha_vec, alpha, part)
%% xfoil and data
load('NACA2412_aeroDB_Abbott_p136.mat');

naca_2412_Re_3p1e6_BL = readtable("naca_2412_Re_3p1e6_BL.txt");
naca_2412_Re_5p7e6_BL = readtable("naca_2412_Re_5p7e6_BL.txt");
naca_2412_Re_8p9e6_BL = readtable("naca_2412_Re_8p9e6_BL.txt");
VariableNames_BL = {'s', 'x','y', 'Ue/vinf','delta_star','theta', 'cf', 'H'};
naca_2412_Re_3p1e6_BL.Properties.VariableNames = VariableNames_BL;
naca_2412_Re_5p7e6_BL.Properties.VariableNames = VariableNames_BL;
naca_2412_Re_8p9e6_BL.Properties.VariableNames = VariableNames_BL;

naca_2412_Re_3p1e6_aseq = readtable("naca_2412_Re_3p1e6_aseq.txt");
naca_2412_Re_5p7e6_aseq = readtable("naca_2412_Re_5p7e6_aseq.txt");
naca_2412_Re_8p9e6_aseq = readtable("naca_2412_Re_8p9e6_aseq.txt");

switch part
    case 'B1 and B2'
        color_vec = generateColorVector(length(Reynolds_list));

        %% plotting

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_upper, Output.(field_name).theta_vec_upper, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_upper, Output.(field_name).theta_vec_upper(Output.(field_name).transition_index_upper*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(1:88)), flip(naca_2412_Re_3p1e6_BL.theta(1:88)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(1:88)), flip(naca_2412_Re_5p7e6_BL.theta(1:88)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(1:88)), flip(naca_2412_Re_8p9e6_BL.theta(1:88)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['$\theta$ for upper surface for $\alpha$ = ' num2str(alpha)], 'Interpreter', 'latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_lower, Output.(field_name).theta_vec_lower, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_lower, Output.(field_name).theta_vec_lower(Output.(field_name).transition_index_lower*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(88:160)), flip(naca_2412_Re_3p1e6_BL.theta(88:160)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(88:160)), flip(naca_2412_Re_5p7e6_BL.theta(88:160)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(88:160)), flip(naca_2412_Re_8p9e6_BL.theta(88:160)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['$\theta$ for lower surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_upper, Output.(field_name).delta_star_upper, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_upper, Output.(field_name).delta_star_upper(Output.(field_name).transition_index_upper*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(1:88)), flip(naca_2412_Re_3p1e6_BL.delta_star(1:88)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(1:88)), flip(naca_2412_Re_5p7e6_BL.delta_star(1:88)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(1:88)), flip(naca_2412_Re_8p9e6_BL.delta_star(1:88)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['$\delta^*$ for upper surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_lower, Output.(field_name).delta_star_lower, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_lower, Output.(field_name).delta_star_lower(Output.(field_name).transition_index_lower*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(88:160)), flip(naca_2412_Re_3p1e6_BL.delta_star(88:160)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(88:160)), flip(naca_2412_Re_5p7e6_BL.delta_star(88:160)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(88:160)), flip(naca_2412_Re_8p9e6_BL.delta_star(88:160)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['$\delta^*$ for lower surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_upper, Output.(field_name).H_vec_upper, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_upper, Output.(field_name).H_vec_upper(Output.(field_name).transition_index_upper*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(1:88)), flip(naca_2412_Re_3p1e6_BL.H(1:88)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(1:88)), flip(naca_2412_Re_5p7e6_BL.H(1:88)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(1:88)), flip(naca_2412_Re_8p9e6_BL.H(1:88)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['H for upper surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_lower, Output.(field_name).H_vec_lower, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_lower, Output.(field_name).H_vec_lower(Output.(field_name).transition_index_lower*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(88:160)), flip(naca_2412_Re_3p1e6_BL.H(88:160)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(88:160)), flip(naca_2412_Re_5p7e6_BL.H(88:160)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(88:160)), flip(naca_2412_Re_8p9e6_BL.H(88:160)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['H for lower surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_upper, Output.(field_name).cf_upper, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_upper, Output.(field_name).cf_upper(Output.(field_name).transition_index_upper*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(1:88)), flip(naca_2412_Re_3p1e6_BL.cf(1:88)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(1:88)), flip(naca_2412_Re_5p7e6_BL.cf(1:88)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(1:88)), flip(naca_2412_Re_8p9e6_BL.cf(1:88)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['$c_f$ for upper surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure;
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            plot(Output.(field_name).x_vec_out_lower, Output.(field_name).cf_lower, 'DisplayName', ['Numerical Re = ', num2str(Re_L)],'Color',color_vec(i,:))
            scatter(Output.(field_name).x_transition_lower, Output.(field_name).cf_lower(Output.(field_name).transition_index_lower*41),'filled','DisplayName', ['Numerical transition Re = ', num2str(Re_L)])
            i = i + 1;
        end
        plot(flip(naca_2412_Re_3p1e6_BL.x(88:160)), flip(naca_2412_Re_3p1e6_BL.cf(88:160)), 'DisplayName', 'Xfoil Re = 3.1*10^6','Color',color_vec(1,:),'LineStyle','--')
        plot(flip(naca_2412_Re_5p7e6_BL.x(88:160)), flip(naca_2412_Re_5p7e6_BL.cf(88:160)), 'DisplayName', 'Xfoil Re = 5.7*10^6','Color',color_vec(2,:),'LineStyle','--')
        plot(flip(naca_2412_Re_8p9e6_BL.x(88:160)), flip(naca_2412_Re_8p9e6_BL.cf(88:160)), 'DisplayName', 'Xfoil Re = 8.9*10^6','Color',color_vec(3,:),'LineStyle','--')
        title(['$c_f$ for lower surface for $\alpha$ = ' num2str(alpha)],'Interpreter','latex')
        xlabel('$\frac{x}{c}$ from stagnation point')
        legend('Location', 'best')
        grid on

        figure
        hold on
        plot(naca_2412_Re_3p1e6_aseq.alpha, naca_2412_Re_3p1e6_aseq.CL, 'r')
        plot(NACA2412_Cl_vs_alpha{1}(:,1),NACA2412_Cl_vs_alpha{1}(:,2),'r--')
        plot(naca_2412_Re_5p7e6_aseq.alpha, naca_2412_Re_5p7e6_aseq.CL, 'b')
        plot(NACA2412_Cl_vs_alpha{2}(:,1),NACA2412_Cl_vs_alpha{2}(:,2),'b--')
        plot(naca_2412_Re_8p9e6_aseq.alpha, naca_2412_Re_8p9e6_aseq.CL, 'k')
        plot(NACA2412_Cl_vs_alpha{3}(:,1),NACA2412_Cl_vs_alpha{3}(:,2),'k--')
        legend('C_L Xfoil Re = 3.1*10^6','C_L experiment Re = 3.1*10^6','C_L Xfoil Re = 5.7*10^6','C_L experiment Re = 5.7*10^6','C_L Xfoil Re = 8.9*10^6','C_L experiment Re = 8.9*10^6','Location','best')
        grid on


        figure
        hold on
        i = 1;
        for Re_L = Reynolds_list
            field_name = ['Re_' num2str(Re_L)];
            scatter(Output.(field_name).CL, Output.(field_name).CD, 'filled','DisplayName',['numerical solution for Re = ',num2str(Re_L) ])
            i = i + 1;
        end
        plot(naca_2412_Re_3p1e6_aseq.CL, naca_2412_Re_3p1e6_aseq.CD, 'r','DisplayName','C_D/C_L Xfoil Re = 3.1*10^6')
        plot(NACA2412_Cd_vs_Cl{1}(:,1),NACA2412_Cd_vs_Cl{1}(:,2),'r--','DisplayName','C_D/C_L experiment Re = 3.1*10^6')
        plot(naca_2412_Re_5p7e6_aseq.CL, naca_2412_Re_5p7e6_aseq.CD, 'b','DisplayName','C_D/C_L Xfoil Re = 5.7*10^6')
        plot(NACA2412_Cd_vs_Cl{2}(:,1),NACA2412_Cd_vs_Cl{2}(:,2),'b--','DisplayName','C_D/C_L experiment Re = 5.7*10^6')
        plot(naca_2412_Re_8p9e6_aseq.CL, naca_2412_Re_8p9e6_aseq.CD, 'k','DisplayName','C_D/C_L Xfoil Re = 8.9*10^6')
        plot(NACA2412_Cd_vs_Cl{3}(:,1),NACA2412_Cd_vs_Cl{3}(:,2),'k--','DisplayName','C_D/C_L experiment Re = 8.9*10^6')
        legend('Location', 'best')
        xlabel('$C_L$','Interpreter','latex')
        ylabel('$C_D$','Interpreter','latex')
        grid on

        figure
        hold on
        plot(naca_2412_Re_3p1e6_aseq.alpha, naca_2412_Re_3p1e6_aseq.CM, 'r')
        plot(NACA2412_Cm025_vs_alpha{1}(:,1),NACA2412_Cm025_vs_alpha{1}(:,2),'r--')
        plot(naca_2412_Re_5p7e6_aseq.alpha, naca_2412_Re_5p7e6_aseq.CM, 'b')
        plot(NACA2412_Cm025_vs_alpha{2}(:,1),NACA2412_Cm025_vs_alpha{2}(:,2),'b--')
        plot(naca_2412_Re_8p9e6_aseq.alpha, naca_2412_Re_8p9e6_aseq.CM, 'k')
        plot(NACA2412_Cm025_vs_alpha{3}(:,1),NACA2412_Cm025_vs_alpha{3}(:,2),'k--')
        legend('C_M Xfoil Re = 3.1*10^6','C_M experiment Re = 3.1*10^6','C_M Xfoil Re = 5.7*10^6','C_M experiment Re = 5.7*10^6','C_M Xfoil Re = 8.9*10^6','C_M experiment Re = 8.9*10^6','Location','best')
        grid on
    case 'B4'
        alpha_vec_Re_3p1e6 = NACA2412_Cl_vs_alpha{1}(:,1);
        alpha_vec_Re_5p7e6 = NACA2412_Cl_vs_alpha{2}(:,1);
        alpha_vec_Re_8p9e6 = NACA2412_Cl_vs_alpha{3}(:,1);
        Cd_alpha_Re_3p1e6 = interp1(NACA2412_Cd_vs_Cl{1}(:, 1), NACA2412_Cd_vs_Cl{1}(:, 2), NACA2412_Cl_vs_alpha{1}(:,2));
        Cd_alpha_Re_5p7e6 = interp1(NACA2412_Cd_vs_Cl{2}(:, 1), NACA2412_Cd_vs_Cl{2}(:, 2), NACA2412_Cl_vs_alpha{2}(:,2));
        Cd_alpha_Re_8p9e6 = interp1(NACA2412_Cd_vs_Cl{3}(:, 1), NACA2412_Cd_vs_Cl{3}(:, 2), NACA2412_Cl_vs_alpha{3}(:,2));

        figure
        hold on
        plot(alpha_vec_Re_3p1e6, Cd_alpha_Re_3p1e6, 'DisplayName', 'Re = 3.1e6 Experiment')
        plot(naca_2412_Re_3p1e6_aseq.alpha, naca_2412_Re_3p1e6_aseq.CD, 'DisplayName', 'Re = 3.1e6 Xfoil')
        for alpha = alpha_vec
            scatter(alpha, Output.Re_3100000.(['Alpha_' num2str(alpha)]).CD,'filled','DisplayName', ['Re = 3.1e6, \alpha = ' num2str(alpha)])
        end
        grid on
        legend show
        title('CD polar for Re  = 3.1e6')


        figure
        hold on
        plot(alpha_vec_Re_5p7e6, Cd_alpha_Re_5p7e6, 'DisplayName', 'Re = 5.7e6 Experiment')
        plot(naca_2412_Re_5p7e6_aseq.alpha, naca_2412_Re_5p7e6_aseq.CD, 'DisplayName', 'Re = 5.7e6 Xfoil')
        for alpha = alpha_vec
            scatter(alpha, Output.Re_5700000.(['Alpha_' num2str(alpha)]).CD,'filled','DisplayName', ['Re = 5.7e6, \alpha = ' num2str(alpha)])
        end
        grid on
        legend show
        title('CD polar for Re  = 5.7e6')

        figure
        hold on
        plot(alpha_vec_Re_8p9e6, Cd_alpha_Re_8p9e6, 'DisplayName', 'Re = 8.9e6 Experiment')
        plot(naca_2412_Re_8p9e6_aseq.alpha, naca_2412_Re_8p9e6_aseq.CD, 'DisplayName', 'Re = 8.9e6 Xfoil')
        for alpha = alpha_vec
            scatter(alpha, Output.Re_8900000.(['Alpha_' num2str(alpha)]).CD,'filled','DisplayName', ['Re = 8.9e6, \alpha = ' num2str(alpha)])
        end
        grid on
        legend show
        title('CD polar for Re  = 8.9e6')


end
