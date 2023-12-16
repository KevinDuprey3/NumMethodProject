function [numDrones,detect] = scouting(vine,t)

    numDrones = 1;
    diam = (0.5 * 20)/10;
    size = pi*(diam/2)^2;
    for i = 1:25
        for j = 1:25

            if vine(i).I(j) >= size
                detect = true;
            else
                detect = false;
            end
        end
    end

end