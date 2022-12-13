function beamFittingLoss = ESMBeamFitter(flightData,geomData,inputChirpData,arrayFullGeomData,desiredSignals,initialBeamData)

%construct the beam pattern from the input
beamPatternData = struct("beamData",[]);
phiSpan = -2*pi:0.1:2*pi;
thetaSpan = -2*pi:0.1:2*pi;

%add the single beam pattern we generate
beamPatternVars = struct("Amp",initialBeamData(1),"Sigx",initialBeamData(2),"Sigy",initialBeamData(3),"Vo",0,"phiSpan",phiSpan,"thetaSpan",thetaSpan);
beamPatternData(1).beamData = beamPatternVars;

%compute the solution for the input beam pattern
fullSolutionStruct = computeESMFieldChirps(flightData,geomData,inputChirpData,beamPatternData,arrayFullGeomData);

%get the specific geometry solution needed
dataSolution = fullSolutionStruct(1).geomSolutions;

%BASIC SQUARED ERROR
%dataSolution = dataSolution(:,1:size(desiredSignals,2));
%meanError = sum(mean(sqrt((desiredSignals - dataSolution).^2)));

% %TOTAL ENVELOPE MEAN SQUARED ERROR
% %for the computed pattern, subtract the actual from desired and compute
% %mean squared error, but do the envelope instead, then compute the sum of
% %the total positive envelope and then compare that to the desired signal
% desiredSigEnvelopes = zeros(size(desiredSignals));
% dataSolutionEnvelopes = zeros(size(dataSolution));
% for i = 1:size(desiredSignals,1)
%     [upd,lod] = envelope(desiredSignals(i,:),30,'peak');
%     desiredSigEnvelopes(i,:) = upd;
%     [ups,los] = envelope(dataSolution(i,:),30,'peak');
%     dataSolutionEnvelopes(i,:) = ups;
% end
% 
% %compute the sums of the envelopes
% sumDes = sum(desiredSigEnvelopes,2);
% sumSol = sum(dataSolutionEnvelopes,2);
% 
% %get the mean error
% meanError = sum(sqrt( (sumDes - sumSol).^2));

%LINED UP ENVELOPE MEAN SQUARED ERROR
%basically do the same envelope stuff as before, but line up the solution
%for each of the thirty microphones, get the flat intervals
%between these is where we can compute the mean mic amplitudes
%get the locations of the solution chirps
diffMicAmps = diff(dataSolution,2,2);

%set the values
oneMicAmps = diffMicAmps;
oneMicAmps(abs(oneMicAmps) > 0) = 1;

%take the diff of this and find where the ones are
onesAtChirps = diff(oneMicAmps,2,2);
onesAtChirps(abs(onesAtChirps) > 0) = 1;

%count the number of ones in each row, they should all be the same
%also get where the chirps are at
solNumChirps = [];
solChirpLocs = [];
for i = 1:size(onesAtChirps,1)
    solNumChirps = [solNumChirps;nnz(onesAtChirps(i,:))];
    %find the chirp locations, if there are no chirps, then set the chirp
    %loation to 1
    chirpLoc = find(onesAtChirps(i,:));
    if length(chirpLoc) > 1
        chirpLoc = chirpLoc(1);
    end
    if (length(chirpLoc) ~= size(solChirpLocs,2)) && ~isempty(solChirpLocs)
        chirpLoc = 1;
    end
    solChirpLocs = [solChirpLocs;chirpLoc];
end

% %get the location of the other chirps
% %for each of the thirty microphones, get the flat intervals
% %between these is where we can compute the mean mic amplitudes
% diffMicAmps1 = diff(desiredSignals);
% 
% %set the values
% oneMicAmps1 = diffMicAmps1;
% oneMicAmps1(abs(oneMicAmps1) > 0) = 1;
% 
% %take the diff of this and find where the ones are
% onesAtChirps1 = diff(oneMicAmps1,2,2);
% onesAtChirps1(abs(onesAtChirps1) > 0) = 1;
% 
% %count the number of ones in each row, they should all be the same
% %also get where the chirps are at
% desiredNumChirps = [];
% desiredChirpLocs = [];
% for i = 1:size(onesAtChirps1,1)
%     desiredNumChirps = [desiredNumChirps;nnz(onesAtChirps1(i,:))];
%     desiredChirpLocs = [desiredChirpLocs;find(onesAtChirps1(i,:))];
% end

%for the solution chirps, shift them all to the right
shiftedSolChirps = zeros(size(dataSolution));
for i = 1:size(dataSolution,1)
    e = size(dataSolution,2) - solChirpLocs(i,1) + 1;
    s1 = dataSolution(i,solChirpLocs(i,1):end);
    shiftedSolChirps(1:size(dataSolution,2) - solChirpLocs(i,1) + 1) = s1;
end

%shiftedSolChirps(1:(size(dataSolution,2) - solChirpLocs(:,1) + 1)) = dataSolution(:,solChirpLocs(:,1):end);

%now concatenate to the desired signal length
shiftedSolChirps = shiftedSolChirps(1:size(desiredSignals,2));

%now get the mean squared error
meanError = sum(mean(sqrt((desiredSignals - shiftedSolChirps).^2)));

%sum the mean squared error and return
beamFittingLoss = meanError;
end

