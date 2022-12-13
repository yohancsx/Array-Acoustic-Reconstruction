function chirps = customChirp(tInput, inputChirp, chirpSpacing)

%set the output array
chirps = zeros(size(tInput));

%create some chirps
st = 0;
for samp = 1:length(tInput)
    st = st + 1;
    if (st == (length(inputChirp) + chirpSpacing)) && (samp + length(inputChirp) < length(chirps))
        chirps(samp:(samp + length(inputChirp) - 1)) = inputChirp(1:end);
        st = 0;
    end
end
end

