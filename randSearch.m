function [NumDrones,ScoutDetect] = randSearch(vine,t,A)
ScoutDetect = 0;
ScoutSpeed = 0.01; %m/s setting to max for now, can use for loop for optimization
NumPlants = ScoutSpeed*3600;
NumDrones = 1;
for j = 1:NumDrones
    ScoutPos = randperm(NumPlants);
    SX=0;
    SY=0;
    ScTime=0;
    for i = 1:NumPlants
        ScDistance = sqrt((SX-vine(ScoutPos(i)).X)^2+(SY-vine(ScoutPos(i)).Y)^2);
        %update time and pos
        SX = vine(ScoutPos(i)).X;
        SY = vine(ScoutPos(i)).Y;
        ScTime = ScTime + ScDistance/ScoutSpeed;

        if ScTime<= 3600 %We can still run the drone
            InfectSize = vine(ScoutPos(i)).L(t)*A;
            ScoutSize = max([(40*ScoutSpeed)^2/4*pi pi]);
            if InfectSize >= ScoutSize
                ScoutDetect = 1; %true; % We found it
                return
            end
        else
            ScoutDetect = 0; %false; % We didnt found it
            break
        end
    end
end
end
