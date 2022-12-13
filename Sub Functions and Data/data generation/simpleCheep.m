function cheeps = simpleCheep(tInput, freq, cheepLen, cheepSpacing)
%create the simple "cheep" function, which simply outputs a set of cheeps
%of a certain length, frequency with spacing determined by the cheep
%spacing (all in terms of samples)

%set the full cheep timeseries
cheeps = zeros(size(tInput));

for s = 1:length(tInput)
    %set up the cheep array
    cheepArray = zeros([1 cheepLen]);
    
    %create the actual cheep
    cheep = sin(2*pi*freq*cheepArray);
    
    %after cheepSpacing samples, then add to the array
    ts = ts + 1;
    if ts == (cheepSpacing + cheepLen)
        cheeps(s:s+cheepLen) = cheep;
        ts = 0;
    end
end

end

