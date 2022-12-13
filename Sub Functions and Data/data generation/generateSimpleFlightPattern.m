function flightPath = generateSimpleFlightPattern(arrayFullGeomData,tSpan,circRad,vecRad,direcType)
%generates a basic three dimensional flight pattern straight through the
%array, also set if we want random directions, somewhat random directions
%or directions that exactly follow the flight path

%set the number of flight path points
numPathPoints = 7;

%generate flight path points
flightPathPoints = [];
for i = 1:numPathPoints + 1
    flightPathPoints = [flightPathPoints;...
        [2*circRad*rand(1,2) - circRad,(i - 1)*(arrayFullGeomData.maxArrayZ/numPathPoints)]];
end

%generate some line direction vectors interpolate between
lineDirectionVectors = [];
for i = 1:numPathPoints + 1
    %create the new vector in the xy plane and between 1-1.2 m long
    newVec = [2*circRad*rand(1,2) - circRad,(arrayFullGeomData.maxArrayZ/numPathPoints)];
    newVec = newVec/norm(newVec);
    newVec = vecRad*rand*newVec;
    lineDirectionVectors = [lineDirectionVectors;newVec];
end


%create the interpolations
interpPoints = linspace(flightPathPoints(1,3),flightPathPoints(end,3),length(tSpan));


%using the interpolations, generate the directions and vectors
%location interpolations
locInterpsX = interp1(flightPathPoints(:,3),flightPathPoints(:,1),interpPoints,'spline');
locInterpsY = interp1(flightPathPoints(:,3),flightPathPoints(:,2),interpPoints,'spline');

%get the vector interpolations
%flight smoothing
vecFlightDirectInterpsX = interp1(flightPathPoints(:,3),lineDirectionVectors(:,1),interpPoints,'spline');
vecFlightDirectInterpsY = interp1(flightPathPoints(:,3),lineDirectionVectors(:,2),interpPoints,'spline');

%add locations together
batLocations = [locInterpsX + vecFlightDirectInterpsX; locInterpsY + vecFlightDirectInterpsY; interpPoints]';
batLocations = batLocations + [arrayFullGeomData.centerPoint(1),arrayFullGeomData.centerPoint(2),0];

%bat direction interpolations
batDirections = [];
%type one, exactly follor the flight path
if direcType == 1
     batDirections = diff(batLocations);
     batDirections = [batDirections;batDirections(end,:)];
%type two, directly folow the flight path
elseif direcType == 2
     batDirections = diff(batLocations);
     batDirections = [batDirections;batDirections(end,:)];
     
     %add a small amount of randomness
     batDirections =  batDirections + [2*0.1*rand(size(batDirections,1),2) - 0.1,zeros([size(batDirections,1),1])];
%type three, completely random
else
    %generate 5 direction vectors where the bat will be pointing at each time
    flightDirecVectors = [];
    for i = 1:numPathPoints + 1
        flightDirecVectors = [flightDirecVectors; rand(1,3)];
    end
    
    %bat pointing direction
    vecDirecInterpsX = interp1(flightPathPoints(:,3),flightDirecVectors(:,1),interpPoints,'spline');
    vecDirecInterpsY = interp1(flightPathPoints(:,3),flightDirecVectors(:,2),interpPoints,'spline');
    vecDirecInterpsZ = interp1(flightPathPoints(:,3),flightDirecVectors(:,3),interpPoints,'spline');
    
    %set the bat direction
    batDirections = [vecDirecInterpsX; vecDirecInterpsY; vecDirecInterpsZ]';
end

%return everything
flightPath = struct('Locations',batLocations,'Directions',batDirections);
end

