function solution = computeESMField(arrayFullGeomData,sourceVars,soundSignal,beamPatternVars,simGeomData,Fs)
%just compute the ESM field for a single time step, in the easiest possible
%way (for now)

%create the sources struct (THIS WE CAN PRECOMPUTE USING THE ARRAY GEOM DATA AND THE SOURCEVARS)
reflectionSourcesStruct = struct("sourceVectors",[],"sourcePoints",[]);
reflectionSourcesStruct = createESMReflectionSources(arrayFullGeomData,sourceVars);
reflectionSourcesStruct.sourceVectors = [reflectionSourcesStruct.sourceVectors; sourceVars.Direction];
reflectionSourcesStruct.sourcePoints = [reflectionSourcesStruct.sourcePoints; sourceVars.Location];

%create the solution struct
solution = struct('solution',[]);

%loop through each geometry
for g = 1:length(simGeomData)
    
    %get the evaluation points
    evalPoints = simGeomData(g).geometryData;
    
    %set the eval points solutions
    evalPointsSolution = [];

    %loop through each point
    for pt = 1:size(evalPoints,1)
        
        %set the single eval point solution
        evalPointSolSources = [];
        
        %loop through each source (FOR NOW JUST DO ONE SOURCE)
        for src = 1:size(reflectionSourcesStruct.sourcePoints,1)
            %NORMALLY THIS WOULD BE TAKEN FROM THE SOURCES STRUCT
            %get the bat orientation and direction 
            batOrientation = reflectionSourcesStruct.sourceVectors(src,:);
            batPos = reflectionSourcesStruct.sourcePoints(src,:);
            singleSourceVars = struct("Location",batPos,"Direction",batOrientation);
            
            %convert the point from cartesian to spherical
            [th,ro,r] = cart2sph(evalPoints(pt,1) - batPos(1),evalPoints(pt,2)  - batPos(2),evalPoints(pt,3) - batPos(3));
            
            %TODO: Modify code for imperfect reflection here (energy loss)
            %for all except the last rouce, which is the initial source

            %generate the beam pattern here
            [beamPatternLin,beamPatternSph] = generateBeamPatternBasic(singleSourceVars,beamPatternVars);
            
            %get the amplitude at that point based on the source
            ampAtPoint = interp2(beamPatternVars.thetaSpan,beamPatternVars.phiSpan,beamPatternLin.amplitude,th,ro);
            
            %add the contribution from the 1/r squared term
            ampAtPoint = ampAtPoint*(1/(r^2));
            
            %calculate the sampling delay
            samplingDelay = r/343;
            
            %calculate the signal delay
            sigDelay = round(samplingDelay*Fs);
            
            %if the signal delay is 0, set it to 1
            if sigDelay == 0
                sigDelay = 1;
            end
            
            %if the delay is too long, then just give the backgorund noise
            %(0) in this case
            if sigDelay > length(soundSignal)
                evalPointSolSources = [evalPointSolSources;0];
            else
                %set the eval point amplitude
                evalPointSolSources = [evalPointSolSources;soundSignal(sigDelay)*ampAtPoint];
            end
        end
        
        %set the eval points solution here
        evalPointsSolution = [evalPointsSolution;sum(evalPointSolSources)];
    end
    
    %add the solution to the evaluation points solution struct
    solution(g).solution = evalPointsSolution;
end

end

