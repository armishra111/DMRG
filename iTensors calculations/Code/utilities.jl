module ITensorCustomFunction
using ITensors
using DelimitedFiles
using ITensors.HDF5

export Expect_N_Values, Expect_NN_Values, Expect_AAdag_Values, Entanglement_Entropy, BoseHubbard_H, DMRG
export Read_PSI_H5, Write_PSI_H5, Make_Directory, Write_To_CSV, Read_From_CSV
export OPEN_BOUND, CLOSE_BOUND

function Expect_N_Values(psi::MPS)
    values = expect(psi, "N")
    return values
end

function Expect_N_Total(psi::MPS)
    values = Expect_N_Values(psi::MPS)
    return sum(values)
end

function Expect_NN_Values(psi::MPS)
    values = correlation_matrix(psi,"N","N")
    return values 
end

function Expect_AAdag_Values(psi::MPS)
    values = correlation_matrix(psi, "A", "Adag")
    return values
end

function Entanglement_Entropy(psi::MPS, b::Int64)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end
    return SvN
end

function BoseHubbard_H(sites, J::Float64, mu::Float64, U::Float64, V::Float64, k::Float64,BC::Int64)
    L = size(sites)[1]
    ampo = OpSum()
    
    for l in 1:(L-1)
        ampo += -J, "A", l, "Adag", l + 1
        ampo += -J, "Adag", l, "A", l + 1
        ampo += V, "N", l, "N", l + 1
    end
    for l in 1:L
        dl = ((L+1)/2-l)
        ampo += -mu, "N", l
        ampo += -U/2, "N", l
        ampo += U/2, "N", l, "N", l
        ampo += k*dl^2, "N", l
    end
    if BC == 1
        ampo += -J, "A", L, "Adag", 1
        ampo += -J, "Adag", L, "A", 1
        ampo += V, "N", L, "N", 1
    end 
    H = MPO(ampo,sites)
    return H
end

function DMRG(H::MPO, sites::Vector{Index{Int64}}, nsweeps::Int64, max_dim, cutoff, noise; psi_guess=nothing, psi_ortho=nothing)
    if isnothing(psi_guess)
        psi_guess = randomMPS(sites,1)
    end
    if isnothing(psi_ortho)
        energy, psi = dmrg(H,psi_guess; nsweeps, max_dim, cutoff, noise)  
    else
        energy, psi = dmrg(H,psi_ortho,psi_guess; nsweeps, max_dim, cutoff, noise)  
    end
    return energy,psi
end

function Read_PSI_H5(filepath)
    f = h5open(filepath,"r")
    psi = read(f,"psi",MPS)
    close(f)
    return psi
end

function Write_PSI_H5(filepath, psi)
    f = h5open(filepath,"w")
    write(f,"psi",psi)
    close(f)
end

function Make_Directory(filepath::String)
    if !isdir(filepath)
        mkdir(filepath)
    else
        println("Directory $filepath already exists.")
    end
end

function Write_To_CSV(filepath, data)
    writedlm(filepath, data , ',')
end

function Read_From_CSV(filepath)
    return readdlm(filepath, ',')
end

OPEN_BOUND = 0
CLOSE_BOUND = 1

end