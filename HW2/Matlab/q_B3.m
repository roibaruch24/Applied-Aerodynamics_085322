%% B3
alpha = 10*pi/180;
A_n_vec = zeros(1,1); 
A_n_vec(1) = alpha -0.0046;
i=1;
while i>0
syms x
exp1 = cos(i*x)*(0.125*cos(x)-0.025);
exp2 = cos(i*x)*(0.0555*cos(x)-0.011);
I1 = int(exp1,x,[0 acos(1/5)]);
I2 = int(exp2,x,[acos(1/5) pi]);
I1 = double(I1);
I2 = double(I2);
A_n_vec(end+1) = (2/pi)*(I1+I2);
if abs(A_n_vec(end)) <= 1*10^-3
    break
end
i=i+1;
end

phi = 10:1:170;
phi = deg2rad(phi);
gamma = zeros(1,length(phi)+1);

for i = 1:length(phi)
    a =phi(i)
    b= rad2deg(phi(i))
    gamma(i) = 2*A_n_vec(1)*(1+cos(phi(i)))/sin(phi(i));
    for j = 1:length(A_n_vec)-1
        gamma(i) = (gamma(i)+2*A_n_vec(j+1)*sin(j*phi(i)));
    end
end

Cp_thinaifoil = 2*gamma;
phi = [phi,pi];
x = 0.5*(1-cos(phi));


figure
plot(x,Cp_thinaifoil)
title('plot of Cp using thin airfoil theory')
xlabel('x/c')
ylabel('\Delta C_p')
grid on



figure
plot(x,Cp_thinaifoil)
hold on
title('comparison of \Delta Cp for thin airfoil and panel method')
xlabel('x/c')
ylabel('\Delta C_p')
grid on
hold off
legend('Thin airfoil')
