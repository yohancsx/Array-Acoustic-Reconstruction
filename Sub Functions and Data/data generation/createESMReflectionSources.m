function reflectionSourcesStruct = createESMReflectionSources(arrayFullGeomData,sourceVars)

%for the array geometry, create the reflections sources struct
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

reflectionSourcesStruct = struct("sourcePoints",reflectedSourcePoints,"sourceVectors",reflectedVectors);

end

