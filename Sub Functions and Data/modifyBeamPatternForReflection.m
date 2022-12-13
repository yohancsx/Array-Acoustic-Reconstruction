function [outputArg1,outputArg2] = modifyBeamPatternForReflection(shouldChange,beamPatternVars,reflectSourceLoc)
%change the spread and total amplitude of the beam pattern for reflection
%if shouldChange is false (0), then we simply return the input as the
%output
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

