function solution = computeESMFieldOptimized(arrayFullGeomData,sourceVars,soundSignal,beamPatternVars,simGeomData,Fs)
%just compute the ESM field for a single time step, in the easiest possible
%way (for now)

%create the sources struct (THIS WE CAN PRECOMPUTE USING THE ARRAY GEOM DATA AND THE SOURCEVARS)
reflectionSourcesStruct = struct("sourceVectors",[],"sourcePoints",[]);
reflectionSourcesStruct = createESMReflectionSources(arrayFullGeomData,sourceVars);
reflectionSourcesStruct.sourceVectors = [reflectionSourcesStruct.sourceVectors; sourceVars.Direction];
reflectionSourcesStruct.sourcePoints = [reflectionSourcesStruct.sourcePoints; sourceVars.Location];

%create the solution struct
solution = struct('solution',[]);

%vectorize the code
%loop through each geometry
for g = 1:length(simGeomData)
    
    %disp("Computing geometry " + string(g) + " solution");
    %get the evaluation points
    evalPoints = simGeomData(g).geometryData;
    
    %set the eval points solutions
    evalPointsSolution = [];
    
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

            %generate the beam pattern here
            [beamPatternSph, beamPatternLin] = generateBeamPatternBasic(singleSourceVars,beamPatternVars);
            
            %get the amplitude at that point based on the source
            %TODO FIX THIS INTERPOLATION, it's not going to interpolate
            %linearly unfortunately
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
            sigDelays(sigDelays > length(soundSignal)) = 1;
            
            %if the delay is too long, then just give the backgorund noise
            %(0) in this case
            singleSol = soundSignal(sigDelays)'.*ampAtPoints;
            %singleSol(isnan(singleSol)) = 0;
            evalPointsSolution = [evalPointsSolution,singleSol];
    end
    %get the final solution
    solution(g).solution = sum(evalPointsSolution,2);
end



