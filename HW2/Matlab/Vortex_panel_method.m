function [Cp_vpm, CL_vpm, CD_vpm, CM_vpm,CL_kutta_vpm, CM_kutta_vpm,  X_mid] = Vortex_panel_method(Input)
     
    CN1 = zeros(Input.N, Input.N);
    CN2 = zeros(Input.N, Input.N);
    CT1 = zeros(Input.N, Input.N);
    CT2 = zeros(Input.N, Input.N);
    A_N = zeros(Input.N+1, Input.N+1);
    A_T = zeros(Input.N, Input.N+1);
    V_tangent = zeros(Input.N,1);
    RHS = zeros(Input.N+1,1);
    Input.alpha =  Input.alpha/57.3;

     X_mid = zeros(1,Input.N);
     Y_mid = zeros(1,Input.N);
     l = zeros(1,Input.N);
     % sin_theta = zeros(1,Input.N);
     % cos_theta = zeros(1,Input.N);

     for i = 1:(length(Input.Y_vec)-1)
        X_mid(i) = (Input.X_vec(i) + Input.X_vec(i+1))/2;
        Y_mid(i) = (Input.Y_vec(i) + Input.Y_vec(i+1))/2;
        l(i) = sqrt((Input.X_vec(i+1)-Input.X_vec(i))^2+(Input.Y_vec(i+1)-Input.Y_vec(i))^2);
        sin_theta(i) = (Input.Y_vec(i+1)-Input.Y_vec(i))/l(i);
        cos_theta(i) = (Input.X_vec(i+1)-Input.X_vec(i))/l(i);
     end
     for i = 1:Input.N
        for j = 1:Input.N
            if i == j
                CN1(i, j) = -1.0;
                CN2(i, j) = 1.0;
                CT1(i, j) = 0.5 * pi;
                CT2(i, j) = 0.5 * pi;
            else
                [sine_value, cosine_value] = sine_cosine_diff (sin_theta(i), sin_theta(j), cos_theta(i), cos_theta(j));
                [sin_i_m_2j,cos_i_m_2j]=sine_cosine_i_m2j_diff(sin_theta(i),sin_theta(j),cos_theta(i),cos_theta(j));
                A = -(X_mid(i) - Input.X_vec(j)) * cos_theta(j) - (Y_mid(i) - Input.Y_vec(j)) * sin_theta(j);
                B = (Input.Y_vec(j)-Y_mid(i))^2+(X_mid(i)-Input.X_vec(j))^2;
                C = sine_value;
                D = cosine_value;
                E = (X_mid(i) - Input.X_vec(j)) * sin_theta(j) - (Y_mid(i) - Input.Y_vec(j)) * cos_theta(j);
                F = log(1.0 + l(j) * (l(j) + 2.0 * A) / B);
                G = atan2(E * l(j), B + A * l(j));
                P = (X_mid(i) - Input.X_vec(j)) * sin_i_m_2j + (Y_mid(i) - Input.Y_vec(j)) * cos_i_m_2j;
                Q = (X_mid(i) - Input.X_vec(j)) * cos_i_m_2j - (Y_mid(i) - Input.Y_vec(j)) * sin_i_m_2j;
                CN2(i, j) = D + 0.5 * Q * F / l(j) - (A * C + D * E) * G / l(j);
                CN1(i, j) = 0.5 * D * F + C * G - CN2(i, j);
                CT2(i, j) = C + 0.5 * P * F / l(j) + (A * D - C * E) * G / l(j);
                CT1(i, j) = 0.5 * C * F - D * G - CT2(i, j);
            end
        end
     end

      for i = 1:Input.N
        A_N(i, 1) = CN1(i, 1);
        A_N(i, Input.N+1) = CN2(i, Input.N);
        A_T(i, 1) = CT1(i, 1);
        A_T(i, Input.N+1) = CT2(i, Input.N);
        for j = 2:Input.N
            A_N(i, j) = CN1(i, j) + CN2(i, j-1);
            A_T(i, j) = CT1(i, j) + CT2(i, j-1);
        end
     end
     A_N(Input.N+1, 1) = 1;
     A_N(Input.N+1, Input.N+1) = 1;
     for j = 2:Input.N
        A_N(Input.N+1, j) = 0;
     end
     for i = 1:Input.N
        [sine_value, ~] = sine_cosine_diff (sin_theta(i), sin(Input.alpha), cos_theta(i), cos(Input.alpha));
        RHS(i) = sine_value;
     end

     gamma = A_N \ RHS;
    
    for i = 1:Input.N
        [~, cosine_value] = sine_cosine_diff (sin_theta(i), sin(Input.alpha), cos_theta(i), cos(Input.alpha));
        for j = 1:Input.N+1
            V_tangent(i) = V_tangent(i) + A_T(i, j) * gamma(j);
        end
        V_tangent(i) = V_tangent(i) + cosine_value;
        Cp_vpm(i) = 1 - V_tangent(i)^2;
    end
Cl=zeros(Input.N,1);
Cd=zeros(Input.N,1);
Cm=zeros(Input.N,1);

for i=1:Input.N
    F(i)=-Cp_vpm(i)*l(i);
    Cl(i)=F(i)*cos_theta(i);
    Cd(i)=F(i)*sin_theta(i);
    Cm(i)=(0.25-X_mid(i))*Cl(i)+(0-Y_mid(i))*Cd(i);
end
CL_vpm=sum(Cl);
CD_vpm=sum(Cd);
CM_vpm=sum(Cm);
 
gamma_Hess_Smith = zeros(Input.N,1);
for i=1:Input.N
    gamma_Hess_Smith(i)=(gamma(i+1)+gamma(i))*l(i)*pi;
end

CL_kutta_temp=(2.*gamma_Hess_Smith)/(Input.V_inf*Input.Chord);
CL_kutta_vpm = sum(CL_kutta_temp);
CM_kutta_temp=sum(CL_kutta_temp.*X_mid')*cos(Input.alpha);
CM_kutta_vpm=(-CM_kutta_temp+CL_kutta_vpm*0.25);

     

end