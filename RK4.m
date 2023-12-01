function x_array = RK4(ode,x_i,t_array)
n = length(x_i);
x_array = zeros(1,n);

%for i=2:n
    h = 1;
    k1 = ode(t_array,x_array(:)); % goes wrong here idk why
    k2 = ode(t_array+0.5*h,x_array(:)+0.5*k1*h);
    k3 = ode(t_array+0.5*h,x_array(:)+0.5*k2*h);
    k4 = ode(t_array+h,x_array(:)+k3*h);
    
    phi = (k1+2*k2+2*k3+k4)/6;
    
    step(:) = phi(:) .* h;
    xi(:) = x_i(:)';
    x_array(:) = xi(:) + step(:);
%end

end