% number of sites
function Bases = Generate_Basis(N)
    Bases = dec2bin(0:2^N-1)' - '0';
    Bases = transpose(Bases);
end

