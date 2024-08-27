function [Output] = HW_3_Part_B(Input)
Input.TE_type = 'blunt'; % can be either blunt or sharp
[Input] = Panel_code_init(Input);
[~, CZ_Hess_Smith, Cx_Hess_Smith, ~,CZ_kutta_Hess_Smith, ~, ~, U_tanget, l] = Hess_Smith(Input);

ue_lower = flip(U_tanget(1:end/2));
ue_upper = U_tanget(end/2+1:end);
sign_change_idx = find(diff(sign(ue_lower)));
ue_upper = [flip(ue_lower(1:sign_change_idx)) ;ue_upper];
ue_lower = ue_lower(sign_change_idx+1:end);

L_lower = flip(l(1:end/2));
L_upper = l(end/2 + 1:end);
L_upper = [flip(L_lower(1:sign_change_idx))';L_upper']';
L_lower = L_lower(sign_change_idx+1:end);

dUe_dx_upper = due_dx(ue_upper, L_upper);
dUe_dx_lower = due_dx(ue_lower, L_lower);

ni_upper = ue_upper/Input.Re_L;
ni_lower = ue_lower/Input.Re_L;

% Upper
Input.L = 0;
Input.x_0 = 1e-6;
Input.theta_0 = sqrt((0.075*ni_upper(1))/(dUe_dx_upper(1)));
Input.du_e_dx = dUe_dx_upper;
Output.theta_vec_upper=[];
Output.x_vec_out_upper=[];
Output.H_vec_upper = [];
Output.delta_star_upper = [];
Output.cf_upper = [];
transition_flag_upper = 0;

for i = 1:length(ue_upper)
    Input.du_e_dx = dUe_dx_upper(i);
    Input.ni = ni_upper(i);
    Input.u_e= ue_upper(i);
    Input.L = Input.L + L_upper(i);
    if transition_flag_upper == 0
        [Output_laminar] = laminar_boundary_layer(Input);
        Output.theta_vec_upper = [Output.theta_vec_upper; Output_laminar.theta_Thwaites];
        Output.x_vec_out_upper = [Output.x_vec_out_upper; Output_laminar.x_l];
        Output.H_vec_upper = [Output.H_vec_upper; Output_laminar.H_Thwaites'];
        Output.delta_star_upper = [Output.delta_star_upper; Output_laminar.delta_star_Thwaites'];
        Output.cf_upper = [Output.cf_upper; Output_laminar.cf_Thwaites'];
        Input.theta_0 = Output_laminar.theta_Thwaites(end);
        Input.x_0 = Input.x_0 + L_upper(i);
        switch Input.transition_type
            case 'free'
                Re_theta = Input.Re_L * Output.theta_vec_upper(end);
                Re_x = Input.Re_L *Output_laminar.x_l(end);
                if Re_theta > 1.174*(1+22400/Re_x)*Re_x^0.46
                    Output.x_transition_upper = Output_laminar.x_l(end);
                    transition_flag_upper = 1;
                    Output.transition_index_upper = i;
                    Input.theta_0 = Output_laminar.theta_Thwaites(end);
                end
            case 'forced'
                if Output.x_vec_out_upper(end) >= Input.transition_point
                    Output.x_transition_upper = Output_laminar.x_l(end);
                    transition_flag_upper = 1;
                    Output.transition_index_upper = i;
                    Input.theta_0 = Output_laminar.theta_Thwaites(end);
                end
        end

    end
    if transition_flag_upper == 1 && Output.transition_index_upper ~= i
        Input.H1_0 = H1_H(Output.H_vec_upper(end));
        [Output_turbulent] = turbulent_boundary_layer(Input);
        Input.theta_0 = Output_turbulent.theta_head(end);
        Input.x_0 = Input.x_0 + L_upper(i);
        Output.theta_vec_upper = [Output.theta_vec_upper; Output_turbulent.theta_head];
        Output.x_vec_out_upper = [Output.x_vec_out_upper; Output_turbulent.x_t];
        Output.H_vec_upper = [Output.H_vec_upper; Output_turbulent.H_head'];
        Output.delta_star_upper = [Output.delta_star_upper; Output_turbulent.delta_star_head'];
        Output.cf_upper = [Output.cf_upper; Output_turbulent.cf_head'];
        if Output.H_vec_upper(end) >= 3
            for j = 1:length(Output.H_vec_upper)
                if Output.H_vec_upper(j) >= 3
                    Output.theta_vec_upper = Output.theta_vec_upper(1:j);
                    Output.x_vec_out_upper = Output.x_vec_out_upper(1:j);
                    Output.x_seperation_upper = Output.x_vec_out_upper(j);
                    Output.H_vec_upper = Output.H_vec_upper(1:j);
                    Output.delta_star_upper = Output.delta_star_upper(1:j);
                    Output.cf_upper = Output.cf_upper(1:j);
                    break
                end
            end
        end
    end
end

% Lower

Input.L = 0;
Input.x_0 = 1e-6;
Input.theta_0 = sqrt((0.075*ni_upper(1))/(dUe_dx_upper(1)));
Input.du_e_dx = dUe_dx_lower;
Output.theta_vec_lower=[];
Output.x_vec_out_lower=[];
Output.H_vec_lower = [];
Output.delta_star_lower = [];
Output.cf_lower = [];
transition_flag_lower = 0;

for i = 1:length(ue_lower)
    Input.du_e_dx = dUe_dx_lower(i);
    Input.ni = ni_lower(i);
    Input.u_e= ue_lower(i);
    Input.L = Input.L + L_lower(i);
    if transition_flag_lower == 0
        [Output_laminar] = laminar_boundary_layer(Input);
        Output.theta_vec_lower = [Output.theta_vec_lower; Output_laminar.theta_Thwaites];
        Output.x_vec_out_lower = [Output.x_vec_out_lower; Output_laminar.x_l];
        Output.H_vec_lower = [Output.H_vec_lower; Output_laminar.H_Thwaites'];
        Output.delta_star_lower = [Output.delta_star_lower; Output_laminar.delta_star_Thwaites'];
        Output.cf_lower = [Output.cf_lower; Output_laminar.cf_Thwaites'];

        Input.theta_0 = Output_laminar.theta_Thwaites(end);
        Input.x_0 = Input.x_0 + L_lower(i);
        Re_theta = Input.Re_L * Output.theta_vec_lower(end);
        Re_x = Input.Re_L *Output_laminar.x_l(end);
        if Re_theta > 1.174*(1+22400/Re_x)*Re_x^0.46
            Output.x_transition_lower = Output_laminar.x_l(end);
            transition_flag_lower = 1;
            Output.transition_index_lower = i;
            Input.theta_0 = Output_laminar.theta_Thwaites(end);
        end
    end
    if transition_flag_lower == 1 && Output.transition_index_lower ~= i
        Input.H1_0 = H1_H(Output.H_vec_lower(end));
        [Output_turbulent] = turbulent_boundary_layer(Input);
        Input.theta_0 = Output_turbulent.theta_head(end);
        Input.x_0 = Input.x_0 + L_lower(i);
        Output.theta_vec_lower = [Output.theta_vec_lower; Output_turbulent.theta_head];
        Output.x_vec_out_lower = [Output.x_vec_out_lower; Output_turbulent.x_t];
        Output.H_vec_lower = [Output.H_vec_lower; Output_turbulent.H_head'];
        Output.delta_star_lower = [Output.delta_star_lower; Output_turbulent.delta_star_head'];
        Output.cf_lower = [Output.cf_lower; Output_turbulent.cf_head'];
        if Output.H_vec_lower(end) >= 3
            for j = 1:length(Output.H_vec_lower)
                if Output.H_vec_lower(j) >= 3
                    Output.theta_vec_lower = Output.theta_vec_lower(1:j);
                    Output.x_vec_out_lower = Output.x_vec_out_lower(1:j);
                    Output.x_seperation_lower = Output.x_vec_out_lower(j);
                    Output.H_vec_lower = Output.H_vec_lower(1:j);
                    Output.delta_star_lower = Output.delta_star_lower(1:j);
                    Output.cf_lower = Output.cf_lower(1:j);
                    break
                end
            end
        end
    end
end
Output.CD = 2*Output.theta_vec_upper(end)*ue_upper(end)^((Output.H_vec_upper(end)+5)/2) + 2*Output.theta_vec_lower(end)*abs(ue_lower(end))^((Output.H_vec_lower(end)+5)/2);
Output.CL = CZ_Hess_Smith/cosd(Input.alpha);
end