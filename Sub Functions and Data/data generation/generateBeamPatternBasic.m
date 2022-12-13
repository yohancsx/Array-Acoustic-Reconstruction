function [beamPatternSph,beamPatternLin] = generateBeamPatternBasic(sourceVars,beamPatternVars)
%generates a basic beam pattern based on the inputs, goes from 0 to 2 pi

%get the source location and direction
sourceLoc = sourceVars.Location;
sourceDirec = sourceVars.Direction;

%get the gaussian beam pattern characteristics
amp = beamPatternVars.Amp;
sigX = beamPatternVars.Sigx;
sigY = beamPatternVars.Sigy;
vo = beamPatternVars.Vo; 
phiSpan = beamPatternVars.phiSpan;
thetaSpan = beamPatternVars.thetaSpan;

%convert the source vector into a basic theta and phi offset
%CHANGED Z AND X THIS WAS THE PROBLEM UR A DUMBASSSS
[th0,phi0,extraR] = cart2sph(sourceDirec(3),sourceDirec(2),sourceDirec(1));

%generate the theta and phi values
[theta,phi] = meshgrid(phiSpan,thetaSpan);

%generate the initial gaussian, offsetting the initial spherical coordinate
%system by the amount specified, 
beamPattern = basicGaus(theta,phi,th0,phi0,amp,sigX,sigY,vo);

%to generate the cartesian, we simply generate an identical beam pattern
%without offset, then rotate all the points manually using rodroguez's
%formula, this makes for a better looking plot
beamPatternBasic = basicGaus(theta,phi,0,0,amp,sigX,sigY,vo);

%get the axis angle rotation matrix
%get the angle between the direction of the bat and the x axis
angleToRotate = acos(dot(sourceDirec,[1;0;0])/(norm(sourceDirec)));

%get the axis to rotate around, cross the x axis with the direction of the
%bat
axisToRotate = cross(sourceDirec,[1;0;0]);
axisToRotate = axisToRotate/norm(axisToRotate);

%get the rotation matrix
rotationMatrix = rotationVectorToMatrix(angleToRotate*axisToRotate);


%get the lobe and convert
[X,Y,Z] = sph2cart(theta,phi,beamPatternBasic);
rotatedLobe = zeros([3 size(X,1) size(X,1)]);
rotatedLobe(1,:,:) = X;
rotatedLobe(2,:,:) = Y;
rotatedLobe(3,:,:) = Z;

%multiply it by the rotation matrix
rotatedLobe = pagemtimes(rotationMatrix,rotatedLobe);

%translate the cartesian beam pattern and return
xBp = squeeze(rotatedLobe(1,:,:)) + sourceLoc(1);
yBp = squeeze(rotatedLobe(2,:,:)) + sourceLoc(2);
zBp = squeeze(rotatedLobe(3,:,:)) + sourceLoc(3);

beamPatternSph = struct("theta",theta,"phi",phi,"amplitude",beamPattern);
beamPatternLin = struct("X",xBp,"Y",yBp,"Z",zBp);
end

