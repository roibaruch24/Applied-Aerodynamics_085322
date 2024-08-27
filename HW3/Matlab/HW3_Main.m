% Applied Aerodynamics HW3 Main script, Roi Baruch
clear all, close all
global Input
%% A
Reynolds_list = [1e6];

for Re_L = Reynolds_list
    Input.Re_L = Re_L;
    HW_3_Part_A(Input)
end
%% B.1 + B.2
clear all
NACA = '2412';
Reynolds_list = [3.1e6 5.7e6 8.9e6];
Input.transition_type = 'free'; % can be either free or forced
alpha_vec = [5]; % Deg
Output = struct();

for Re_L = Reynolds_list
    for alpha = alpha_vec
        Input.Re_L = Re_L;
        Input.alpha = alpha;
        Input.NACA = NACA;
        Re_field_name = ['Re_' num2str(Re_L)];
        Output.(Re_field_name) = HW_3_Part_B(Input);
    end
end
%% plot B.1 + B.2
part = 'B1 and B2';
for alpha = alpha_vec
    Input.alpha = alpha;
    Input.NACA = '2412';
    HW_3_Part_B_Plotter(Output, Reynolds_list,alpha_vec, alpha, part);
end
%% B3
clear all, clc;
NACA = '2412';
Reynolds_list = [3.1e6];
Input.transition_type = 'forced'; % can be either free or forced
Input.transition_point = 0.01;    % forcd transition point
alpha_vec = [1:1:10 12:20];             % Deg
Output = struct();

for Re_L = Reynolds_list
    for alpha = alpha_vec
        Input.Re_L = Re_L;
        Input.alpha = alpha;
        Input.NACA = NACA;
        % Re_field_name = ['Re_' num2str(alpha)];
        Alpha_field_name = ['Alpha_' num2str(alpha)];
        Output.(Alpha_field_name) = HW_3_Part_B(Input);
        disp(['alpha = ', num2str(alpha), ', CL = ', num2str(Output.(Alpha_field_name).CL), ', X separation  = ', num2str(Output.(Alpha_field_name).x_seperation_upper)])

    end
end
%% B4
clear all, clc;
part = 'B4';
Reynolds_list = [3.1e6 5.7e6 8.9e6];
NACA = '2412';
Input.transition_type = 'free'; % can be either free or forced
alpha_vec = [1:5];             % Deg
Output = struct();

for Re_L = Reynolds_list
    for alpha = alpha_vec
        Input.Re_L = Re_L;
        Input.alpha = alpha;
        Input.NACA = NACA;
        Re_field_name = ['Re_' num2str(Re_L)];
        Alpha_field_name = ['Alpha_' num2str(alpha)];
        Output.(Re_field_name).(Alpha_field_name) = HW_3_Part_B(Input);
    end
end

Input.NACA = '2412';
HW_3_Part_B_Plotter(Output, Reynolds_list, alpha_vec, alpha, part);
