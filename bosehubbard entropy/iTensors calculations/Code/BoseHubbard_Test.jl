include("utilities.jl")
using .ITensorCustomFunction
using ITensors 

function run()
    L=4
    J,mu,U,V,k,BC = 0.5, 0.5, 1.0, 0.0, 0.0, 0
    sites = siteinds("Boson",L; conserve_qns=false, dim=6)
    psi0 = randomMPS(sites,1)

    H = BoseHubbard_H(sites,J, mu,U, V, k, BC)

    nsweeps = 60
    maxdim = [10,20,30,40,50,80,100]
    cutoff = [1E-13]
    noise = [1E-11]
    @time energy0, psi0 = dmrg(H,psi0; nsweeps, maxdim, cutoff, noise)
end
run()