function [t,y] = runge_kutta_four(odefun,t,y0)

phi = zeros(length(y0),1); % initialize derivative vector

k1 = zeros(length(y0),length(t));
k2 = zeros(length(y0),length(t));
k3 = zeros(length(y0),length(t));
k4 = zeros(length(y0),length(t));
tempParams = params;
y = zeros(length(y0),length(t));
initDay = day;

y(1,1) = y0(1);
y(2,1) = y0(2);
y(3,1) = y0(3);
y(4,1) = y0(4);
y(5,1) = y0(5);


for j = 1 : length(y0)
        
    for i = 1 : length(tspan)
    
        h = t(i+1) - t(i);
        
        day = t(i);
        k1(j,i) = odefun{i}(y(:,i));
        tempParams(j) = y(j,i) + 0.5 * h * k1;
        day = initDay + 0.5 * h;
        k2(j,i) = odefun{i}(tempParams);
        tempParams(j) = y(j,i) + 0.5*k2*h;
        day = initDay + 0.5 * h;
        k3(j,i) = odefun{i}(tempParams);
        tempParams(j) = y(j,i) + k3 * h;
        day = initDay + h;
        k4(j,i) = odefun{i}(tempParams);

    
        phi(j) = (h/6) * (k1(j,i) + 2 * k2(j,i) + 2 * k3(j,i) + k4(j,i));

        y(j,i+1) = y(j,i) + phi(j) * h;

    end


end

end