function [Output] = turbulent_boundary_layer(Input)

% Solve the ODE
[Output.x_t, p] = ode45(@head_ode, [Input.x_0 Input.L], [Input.H1_0, Input.theta_0]);

Re_x_t = zeros(1, length(Output.x_t));
H_head = zeros(1, length(Output.x_t));
cf_head = zeros(1, length(Output.x_t));
delta_star_head = zeros(1, length(Output.x_t));
delta_star_Schlichting = zeros(1, length(Output.x_t));
theta_Schlichting = zeros(1, length(Output.x_t));
cf_Schlichting = zeros(1, length(Output.x_t));
H_Schlichting = zeros(1, length(Output.x_t));
theta_head = p(:, 2);

% Calculate all the parameters
for i = 1:length(Output.x_t)
    Re_x_t(i) = (Input.u_e * Output.x_t(i)) / Input.ni;
    H_head(i) = H_H1(p(i, 1));
    cf_head(i) = (0.246 * 10^(-0.678 * H_head(i))) * ((Input.u_e * theta_head(i)) / Input.ni)^(-0.268);
    delta_star_head(i) = H_head(i) * theta_head(i);
    delta_star_Schlichting(i) = 0.046 * Output.x_t(i) / (Re_x_t(i)^0.2);
    theta_Schlichting(i) = 0.036 * Output.x_t(i) / (Re_x_t(i)^0.2);
    cf_Schlichting(i) = 0.0592 / (Re_x_t(i)^0.2);
    H_Schlichting(i) = delta_star_Schlichting(i)/theta_Schlichting(i);
end

Cf_total_head = trapz(Output.x_t, cf_head)/Input.L;
Cf_total_Schlichting = 0.074/Input.Re_L^(0.2);

Output.Re_x_t = Re_x_t;
Output.H_head = H_head;
Output.cf_head = cf_head;
Output.delta_star_head = delta_star_head;
Output.delta_star_Schlichting = delta_star_Schlichting;
Output.theta_Schlichting = theta_Schlichting;
Output.cf_Schlichting = cf_Schlichting;
Output.theta_head = theta_head;
Output.H_Schlichting = H_Schlichting;
Output.Cf_total_Schlichting = Cf_total_Schlichting;
Output.Cf_total_head = Cf_total_head;
    function head_ode = head_ode(~,p)
            H=H_H1(p(1));
            cf = (0.246*10^(-0.678*H))*(Input.u_e*p(2)/Input.ni)^(-0.268);
            head_ode = [(0.0306/p(2))*(p(1)-3)^-0.6169-p(1)/p(2)*(0.5*cf-p(2)/Input.u_e*(H+2)*Input.du_e_dx)-p(1)/Input.u_e*Input.du_e_dx;
                         0.5*cf-p(2)/Input.u_e*(H+2)*Input.du_e_dx];
            
    end
    function H =H_H1(H1)
        if H1 <= 3.3
            H = 3;
        elseif 3.3 < H1 && H1 < 5.3
            H = 1.1538*(H1 -3.3)^-0.326 + 0.6778;
        else
            H = 0.86*(H1-3.3)^-0.777 + 1.1;
        end
    end
end