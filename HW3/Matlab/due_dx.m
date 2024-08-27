function dUe_dx = due_dx(u_vec, x_vec)

    dUe_dx = zeros(1,length(u_vec));   
    for i = 2:length(x_vec)-1
        dUe_dx(i) = (u_vec(i+1) - u_vec(i-1)) / (0.5*x_vec(i+1)+ x_vec(i)+0.5*x_vec(i-1));
    end
    dUe_dx(1) = (u_vec(2) - u_vec(1)) / (0.5*x_vec(2) + 0.5*x_vec(1));
    dUe_dx(end) = (u_vec(end) - u_vec(end-1)) / (0.5*x_vec(end) +0.5*x_vec(end-1));
end