function [Cp_Hess_Smith, CL_Hess_Smith, CD_Hess_Smith, CM_Hess_Smith,CL_kutta_Hess_Smith, CM_kutta_Hess_Smith, X_mid, U_tanget,l] = Hess_Smith(Input)
    

    X_mid = zeros(1,Input.N);
    Y_mid = zeros(1,Input.N);
    l = zeros(1,Input.N);
    sin_theta = zeros(1,Input.N);
    cos_theta = zeros(1,Input.N);
    A_mat = zeros(Input.N+1,Input.N+1);
    B_vec = zeros(1,Input.N+1);
    beta_ij = zeros(Input.N,Input.N);
    Input.alpha =  Input.alpha/57.3;
    
    for i = 1:(length(Input.Y_vec)-1)
        X_mid(i) = (Input.X_vec(i) + Input.X_vec(i+1))/2;
        Y_mid(i) = (Input.Y_vec(i) + Input.Y_vec(i+1))/2;
        l(i) = sqrt((Input.X_vec(i+1)-Input.X_vec(i))^2+(Input.Y_vec(i+1)-Input.Y_vec(i))^2);
        sin_theta(i) = (Input.Y_vec(i+1)-Input.Y_vec(i))/l(i);
        cos_theta(i) = (Input.X_vec(i+1)-Input.X_vec(i))/l(i);
    end
    
    for i = 1:Input.N
        for j = 1:Input.N
            if i==j
                beta_ij(i,j) = pi;
                r_ij = sqrt((X_mid(i)-Input.X_vec(j))^2+(Y_mid(i)-Input.Y_vec(j))^2);
                r_ij_p1 = sqrt((X_mid(i)-Input.X_vec(j+1))^2+(Y_mid(i)-Input.Y_vec(j+1))^2);
    
            else
                r_ij = sqrt((X_mid(i)-Input.X_vec(j))^2+(Y_mid(i)-Input.Y_vec(j))^2);
                r_ij_p1 = sqrt((X_mid(i)-Input.X_vec(j+1))^2+(Y_mid(i)-Input.Y_vec(j+1))^2);
                
                beta_num = (Input.X_vec(j)-X_mid(i))*(Input.Y_vec(j+1)-Y_mid(i))-(Input.Y_vec(j)-Y_mid(i))*(Input.X_vec(j+1)-X_mid(i));
                
                
                beta_denum = (Input.X_vec(j)-X_mid(i))*(Input.X_vec(j+1)-X_mid(i))+(Input.Y_vec(j)-Y_mid(i))*(Input.Y_vec(j+1)-Y_mid(i));
    
                beta_ij(i,j) = atan2(beta_num,beta_denum);
            end
                % A mat up to Input.N*Input.N
                [sine_value, cosine_value] = sine_cosine_diff (sin_theta(i), sin_theta(j), cos_theta(i), cos_theta(j));
                A_mat(i,j) = (sine_value*log(r_ij_p1/r_ij) + cosine_value*beta_ij(i,j))/(2*pi);
                
                % A mat i, Input.N+1
                A_mat(i,Input.N+1) = A_mat(i,Input.N+1) + (log(r_ij_p1/r_ij)*cosine_value - beta_ij(i,j)*sine_value)/(2*pi);
               
        end
    end
    
    for j=1:Input.N   
       for k=[1,Input.N]
           [sine_value, cosine_value] = sine_cosine_diff(sin_theta(k),sin_theta(j),cos_theta(k),cos_theta(j));
    
           r_kj=sqrt((Y_mid(k)-Input.Y_vec(j))^2+(X_mid(k)-Input.X_vec(j))^2);
           r_Kj_p1=sqrt((Y_mid(k)-Input.Y_vec(j+1))^2+(X_mid(k)-Input.X_vec(j+1))^2);
           
           A_mat(Input.N+1,j) = A_mat(Input.N+1,j) + (beta_ij(k,j)*sine_value-log(r_Kj_p1/r_kj)*cosine_value)/(2*pi);%i,Input.N+1 argument
           A_mat(Input.N+1,Input.N+1) = A_mat(Input.N+1,Input.N+1) + (((log(r_Kj_p1/r_kj)*sine_value)+beta_ij(k,j)*cosine_value))/(2*pi); %Input.N+1,Input.N+1 argument
       end
    end
    
    % B vec 1:Input.N
    for i =1:Input.N
        [sine_value, ~] = sine_cosine_diff (sin_theta(i), sin(Input.alpha), cos_theta(i), cos(Input.alpha));
        B_vec(i) = Input.V_inf*sine_value;
    end
    
    % B vec Input.N+1
    [~, cosine_value_1] = sine_cosine_diff (sin_theta(1), sin(Input.alpha), cos_theta(1), cos(Input.alpha));
    [~, cosine_value_N] = sine_cosine_diff (sin_theta(Input.N), sin(Input.alpha), cos_theta(Input.N), cos(Input.alpha));
    B_vec(Input.N+1) = -Input.V_inf*cosine_value_1-Input.V_inf*cosine_value_N;
    
    A_mat=inv(A_mat);
    q_vec=A_mat*B_vec';
    % q_vec = A_mat\B_vec;
    
    U_tanget = zeros(Input.N,1);
    for  i = 1:Input.N
        for  j = 1:Input.N
            [sine_value, cosine_value] = sine_cosine_diff (sin_theta(i), sin_theta(j), cos_theta(i), cos_theta(j));
            
            r_ij=sqrt((Y_mid(i)-Input.Y_vec(j))^2+(X_mid(i)-Input.X_vec(j))^2);
            r_ij_p1=sqrt((Y_mid(i)-Input.Y_vec(j+1))^2+(X_mid(i)-Input.X_vec(j+1))^2);
    
            U_tanget(i)=U_tanget(i)+(q_vec(j)/(2*pi))*((sine_value)*beta_ij(i,j)-cosine_value*log(r_ij_p1/r_ij))+...
            (q_vec(Input.N+1)/(2*pi))*(sine_value*log(r_ij_p1/r_ij)+cosine_value*beta_ij(i,j));
        end
        [~,cosine_value]=sine_cosine_diff(sin_theta(i),sin(Input.alpha),cos_theta(i),cos(Input.alpha));
        U_tanget(i)=U_tanget(i)+Input.V_inf*cosine_value;
    end
    Cp_Hess_Smith = 1-(U_tanget.^2)/(Input.V_inf^2);
    Cl=zeros(Input.N,1);
    Cd=zeros(Input.N,1);
    Cm=zeros(Input.N,1);
    F=zeros(Input.N,1);
    
    
    for i=1:Input.N
    
        F(i)=-Cp_Hess_Smith(i)*l(i);
        Cl(i)=F(i)*cos_theta(i);
        Cd(i)=F(i)*sin_theta(i);
        Cm(i)=(0.25-X_mid(i))*Cl(i)+(0-Y_mid(i))*Cd(i);
        
    end
    CL_Hess_Smith=sum(Cl);
    CD_Hess_Smith=sum(Cd);
    CM_Hess_Smith=sum(Cm);
     
    gamma=l.*q_vec(end);
    CL_kutta_temp=(2*gamma)/(Input.V_inf*Input.Chord);
    CL_kutta_Hess_Smith = sum(CL_kutta_temp);
    CM_kutta_temp=sum(CL_kutta_temp.*X_mid)*cos(Input.alpha);
    CM_kutta_Hess_Smith=(-CM_kutta_temp+CL_kutta_Hess_Smith*0.25)/pi;
end
