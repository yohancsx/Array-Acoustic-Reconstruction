function simGeomOutput = computeBasicESMFieldSingleStep(arrayFullGeomData,sourceVars,soundSignal,beamPatternVars,simGeomData,Fs)
%Compute the acoustic field ESM for a single timestep, takes the input as
%the array geometry (internal walls not supported), source variables
%(location and direction), the acoustic signal, the beam pattern
%(spherical form) and the geometry data for which to simulate for

%set up the main/core source
sources = struct("sourceLocation",[],"sourceVector",[],"beamPattern",[],"sourceSignal",[]);
sources(1).sourceLocation = sourceVars.Location;
sources(1).sourceVector = sourceVars.Direction;
sources(1).sourceSignal = soundSignal;
%generate beam pattern here
[beamPatternSph, beamPatternLin] = generateBeamPatternBasic(sourceVars,beamPatternVars);
sources(1).beamPattern = beamPatternSph;

%TODO: Create function for this
%first compute the reflected sources and their beam patterns
%compute the positions of the reflected sources
reflectedSourcePoints = [];
for i = 1:size(arrayFullGeomData.arrayWallNorms,1)
    reflectedPoint = dot(sourceVars.Location - arrayFullGeomData.rectangleCenters(i,:),arrayFullGeomData.arrayWallNorms(i,:))*(-2*arrayFullGeomData.arrayWallNorms(i,:)) + sourceVars.Location;
    reflectedSourcePoints = [reflectedSourcePoints; reflectedPoint];
end

%compute the reflection vectors of the sources
reflectedVectors = [];
for i = 1:length(reflectedSourcePoints)
    reflecVec = sourceVars.Direction - 2*dot(sourceVars.Direction,arrayFullGeomData.arrayWallNorms(i,:)).*arrayFullGeomData.arrayWallNorms(i,:);
    reflectedVectors = [reflectedVectors;reflecVec];
end

%add to the source vector struct, making sure to recompute the beam pattern
%direction 
for i = 2:length(reflectedSourcePoints)
    sources(i).sourceLocation = reflectedSourcePoints(i,:);
    sources(i).sourceVector = reflectedVectors(i,:);
    newBeamSourceVars = struct("Location",reflectedSourcePoints(i,:),"Direction",reflectedVectors(i,:));
    %generate beam pattern here (note, this is where we could edit to add
    %in some interesting wall charcteristics)
    [beamPatternSph, beamPatternLin] = generateBeamPatternBasic(newBeamSourceVars,beamPatternVars);
    sources(i).beamPattern = beamPatternSph;
    sources(i).sourceSignal = soundSignal;
end


%for each of the sources, compute the addition to the field from each
%source

%struct to hold the source data
preCompSourceData = struct("sourceDistances",[],"ampDecay",[],"phaseDelay",[],"beamPattDecay",[]);

%loop through each of the sources
for src = 1:size(sources,1)
    
    %for each of the solution geometries
    for geom = 1:length(simGeomData)
        
        %get the array geometry from the struct
        arrayGeomPoints = simGeomData(geom).geometryData;
        
        %translate the points such that the source location is at zero
        arrayGeomPoints = arrayGeomPoints - sources(src).sourceLocation;
    
        %now we convert from the rotated points into spherical 
        [arrPtsTh,arrPtsRo,arrPtsR] = cart2sph(arrayGeomPoints(:,1),arrayGeomPoints(:,2),arrayGeomPoints(:,3));

        %now we compute the beam pattern decay for each theta and ro value
        %using an interpolation
        beamPattDecay = interp2(beamPatternVars.thetaSpan,beamPatternVars.phiSpan,1.5*sources(src).beamPattern.amplitude,arrPtsTh,arrPtsRo);

        %compute the distance of each point to the source
        sourceDistances = vecnorm(arrayGeomPoints - sources(src).sourceLocation,2,2);

        %compute the radial amplitude decay for each point
        pointAmplitudeDecay = 1./(sourceDistances.^2);

        %compute the phase delay of each of the points based on the
        %distance, this is nothing but the distance divided by the speed of
        %sound to find the time delay, then divide by the sampling rate to
        %find the sample delay, note that each point we start at that
        %negative index for the sample delay, then increase the sample for
        %each time step
        phaseDelay = round(-1*(sourceDistances/331.5)/(1/Fs));

        %now save all the data
        preCompSourceData(src,geom).sourceDistances = sourceDistances;
        preCompSourceData(src,geom).ampDecay = pointAmplitudeDecay;
        preCompSourceData(src,geom).phaseDelay = phaseDelay;
        preCompSourceData(src,geom).beamPattDecay = beamPattDecay;
    end
end

%create the output struct
simGeomOutput = struct("solution",[]);

%now actually add the contribution of each source to the acoustic field
%for each of the timesteps
for geom = 1:length(simGeomData)
    
    %set up the solution points vector
    solPointsVector = zeros([size(simGeomData(geom).geometryData,1),1]);
    
    %for each of the sources
    for src = 1:size(sources,1)
        
        %TODO: SOMETHING IS PROBABLY WRONG WITH THE PHASE DELAY STUFF
        %compute the signal based on the phase delay from that source
        pointDelays = preCompSourceData(src,geom).phaseDelay;
        pointDelays(pointDelays <= 0) = 1;
        pointDelays(pointDelays > length(sources(src).sourceSignal)) = 1;
        
        %the signal value at each location
        sigValsAtLocation = sources(src).sourceSignal(pointDelays);
        
        %compute the contrubution from each source 
        sourceContrib = preCompSourceData(src,geom).ampDecay.*sigValsAtLocation'.*preCompSourceData(src,geom).beamPattDecay;
        
        %add to the total array
        solPointsVector = solPointsVector + sourceContrib;
        
        
    end
   
    %add to the solution points struct
    simGeomOutput(geom).solution = solPointsVector;

end

end

