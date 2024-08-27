function [Output] = laminar_boundary_layer(Input)
[x_l, theta_Thwaites] = ode45(@theta_ode, [Input.x_0 Input.L], Input.theta_0);
% calculate all the parameters

Re_x_l = zeros(1,length(x_l));
theta_Blasius = zeros(1,length(x_l));
cf_Blasius = zeros(1,length(x_l));
delta_star_Blasius = zeros(1,length(x_l));
delta_star_Thwaites = zeros(1,length(x_l));
cf_Thwaites = zeros(1,length(x_l));
lambda_Thwaites = zeros(1, length(x_l));
l = zeros(1, length(x_l));
H_Thwaites = zeros(1, length(x_l));
transition_flag = zeros(1, length(x_l));

for i = 1:length(x_l)
    Re_x_l(i) = (Input.u_e*x_l(i))/Input.ni;

    theta_Blasius(i) = (0.664*x_l(i))/sqrt(Re_x_l(i));
    cf_Blasius(i) = 0.664/sqrt(Re_x_l(i));
    delta_star_Blasius(i) = 1.7208*x_l(i)/sqrt(Re_x_l(i));
    
    lambda_Thwaites(i) = Input.du_e_dx*theta_Thwaites(i)^2/Input.ni;
    if lambda_Thwaites(i) > 0.1
        lambda_Thwaites(i) = 0.1;
        l(i) = 0.22 + 1.57*lambda_Thwaites(i)-1.8*lambda_Thwaites(i)^2;
        H_Thwaites(i) = 2.61 - 3.75*lambda_Thwaites(i)-5.24*lambda_Thwaites(i)^2;
    elseif lambda_Thwaites(i)<= 0.1 && lambda_Thwaites(i) >= 0 
        l(i) = 0.22 + 1.57*lambda_Thwaites(i)-1.8*lambda_Thwaites(i)^2;
        H_Thwaites(i) = 2.61 - 3.75*lambda_Thwaites(i)-5.24*lambda_Thwaites(i)^2;
    elseif lambda_Thwaites(i) < 0 && lambda_Thwaites(i) >-0.1
        l(i) = 0.22 + 1.402*lambda_Thwaites(i) + (0.018*lambda_Thwaites(i)/(0.107+lambda_Thwaites(i)));
        H_Thwaites(i) = 2.088 + 0.0731/(0.14+lambda_Thwaites(i));
    else 
        H_Thwaites(i) = 0;
        l(i) = 0;
        transition_flag(i) = 1;
    end
    delta_star_Thwaites(i) = H_Thwaites(i)*theta_Thwaites(i);
    cf_Thwaites(i) = ((2*Input.ni*l(i))/Input.u_e)/theta_Thwaites(i);    
end

Cf_total_Blasius = 1.328/sqrt(Input.Re_L);
Cf_total_Thwaites = trapz(x_l, cf_Thwaites)/Input.L;
H_Blasius = delta_star_Blasius./theta_Blasius;

Output.theta_Thwaites = theta_Thwaites;
Output.x_l = x_l;
Output.theta_Blasius = theta_Blasius;
Output.cf_Blasius = cf_Blasius;
Output.delta_star_Blasius = delta_star_Blasius;
Output.delta_star_Thwaites = delta_star_Thwaites;
Output.cf_Thwaites = cf_Thwaites;
Output.lambda_Thwaites = lambda_Thwaites;
Output.l = l;
Output.H_Thwaites = H_Thwaites;
Output.transition_flag = transition_flag;
Output.Cf_total_Blasius = Cf_total_Blasius;
Output.Cf_total_Thwaites = Cf_total_Thwaites;
Output.H_Blasius = H_Blasius;


function dtheta_dx = theta_ode(~, theta)
    dtheta_dx = (0.225 * Input.ni) / (Input.u_e * theta) - 3 * theta * Input.du_e_dx / Input.u_e;
end

end