function fullSolutionStruct = computeESMFieldChirps(flightData,geomData,inputChirpData,beamPatternData,arrayFullGeomData)
%Computes the ESM field for a series of chirps

%parse the needed data for the function from the input structs
%chirp data
simTimeSteps = inputChirpData.simTimesteps;
Fs = inputChirpData.Fs;
chirpSpacing = inputChirpData.chirpSpacing;
chirpLen = inputChirpData.chirpLen;
totNumChirps = inputChirpData.totNumChirps;
outputSignalOverTime = inputChirpData.outputSignalOverTime;

%flight data
batPosOverTime = flightData.Locations;
batRotationOverTime = flightData.Directions;


%create the full solution struct
fullSolutionStruct = struct("geomSolutions",[],"geometryData",[]);

%create the solution time arrays, which have one large time array for each
%data point in the solution
for g = 1:length(geomData)
    
    %set the data in the solution
    solData = zeros(size(geomData(g).geometryData,1),round(length(simTimeSteps)*1.2));
    
    %set the geometry solution
    fullSolutionStruct(g).geomSolutions = solData;
    
    %set the geometry for each 
    fullSolutionStruct(g).geometryData = geomData(g).geometryData;
    
end

%in this case, we change the second geometry point at each timestep, to a
%point which follows the bat position, initially just set this

%save the t values
tValues = [];

%for each timestep, compute the solution struct
for chirpStep = 1:(totNumChirps) 
    %change this eventually, rn we can have an out of bounds exception
    
    %get the index value
    t = floor((chirpSpacing + chirpLen)*Fs*(chirpStep - 1) + 1);
    tValues = [tValues;t];
    
    %get the chirp at that time
    chirpAtTime = outputSignalOverTime(t:t + (chirpLen)*Fs);
    
    %generate the source vars at this time
    sourceVars = struct("Location",batPosOverTime(t,:),"Direction",batRotationOverTime(t,:));
    
    %get the beam pattern vars at this time
    beamPatternVars = beamPatternData(chirpStep).beamData;
    
    %compute the solution    
    
    %generate the reflection location and directions for the timestep   
    reflectionSourcesStruct = struct("sourceVectors",[],"sourcePoints",[]);
    reflectionSourcesStruct = createESMReflectionSources(arrayFullGeomData,sourceVars);
    reflectionSourcesStruct.sourceVectors = [reflectionSourcesStruct.sourceVectors; sourceVars.Direction];
    reflectionSourcesStruct.sourcePoints = [reflectionSourcesStruct.sourcePoints; sourceVars.Location];

    %for each geometry
    for g = 1:length(geomData)
        
        %get the geometry
        evalPoints = geomData(g).geometryData;
        
        %in this case, for the second geometry, we can change the position
        %over time, also store this in the geometry data
        if g == 2
            evalPoints = batPosOverTime(t,:) + [0 0.1 0];
            
            fullSolutionStruct(2).geometryData = [fullSolutionStruct(2).geometryData;evalPoints];
        end
        
        %for each source
        for src = 1:size(reflectionSourcesStruct.sourcePoints,1)
            
            %get the bat orientation and direction 
            batOrientation = reflectionSourcesStruct.sourceVectors(src,:);
            sourcePos = reflectionSourcesStruct.sourcePoints(src,:);
            singleSourceVars = struct("Location",sourcePos,"Direction",batOrientation);
            
            %convert the point from cartesian to spherical
            %ALSO FIX THIS INTERPOLATION
            [th,elevation,r] = cart2sph((evalPoints(:,3) - sourcePos(3)),(evalPoints(:,2)  - sourcePos(2)),(evalPoints(:,1) - sourcePos(1)));
            
            %TODO: Modify code for imperfect reflection here (energy loss)
            %for all except the last rouce, which is the initial source

            %generate the beam pattern here, or take the vars from
            %pre-generated data
            [beamPatternSph, beamPatternLin] = generateBeamPatternBasic(singleSourceVars,beamPatternVars);
            
            %get the amplitude at that point based on the source
            %ampAtPoints = interp2(beamPatternVars.thetaSpan,beamPatternVars.phiSpan,beamPatternLin.amplitude,th,ro);
            ampAtPoints = interp2(beamPatternSph.theta,beamPatternSph.phi,beamPatternSph.amplitude,th,elevation);
            ampAtPoints(isnan(ampAtPoints)) = 0;
            
            %add the contribution from the 1/r squared term
            r(r == 0) = 0.0001;
            ampAtPoints = ampAtPoints.*(1./(r.^2));
            
            %calculate the sampling delay
            samplingDelays = r./343;
            
            %calculate the signal delay
            sigDelays = round(samplingDelays*Fs);
            
            %if the signal delay is 0, set it to 1
            sigDelays(sigDelays == 0) = 1;
            sigDelays(sigDelays > length(outputSignalOverTime)) = 1;
            
            %for each of the geom points
            for point = 1:length(ampAtPoints)
                
                %get the single time sampling delay
                singleSampDelay = sigDelays(point);
                
                %set the single time amplitude modifier
                ampMultiplier = ampAtPoints(point);
                
                %at the correct location in each point's solution time array 
                %add the portion of the input (after amplitude shifting)
                fullSolutionStruct(g).geomSolutions(point,(t+singleSampDelay):(t+singleSampDelay)+chirpLen*Fs) = ...
                    fullSolutionStruct(g).geomSolutions(point,(t+singleSampDelay):(t+singleSampDelay)+chirpLen*Fs) + chirpAtTime*ampMultiplier;
                
            end %end geometry point for
            
        end %end source for   
        
    end %end geometry for
end
end

