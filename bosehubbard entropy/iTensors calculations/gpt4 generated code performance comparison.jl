using ITensors

function bose_hubbard_dmrg(N::Int, t::Float64, U::Float64, μ::Float64, max_bosons::Int, max_sweeps::Int, m::Int)
  # Initialize the site set with Bose-Hubbard site types
  sites = siteinds("BoseHubbard", N; conserve_Nf=true, max_bosons=max_bosons)

  # Create the Bose-Hubbard Hamiltonian MPO
  ampo = AutoMPO()
  for j in 1:N-1
    ampo += -t, "Bdag", j, "B", j + 1
    ampo += -t, "Bdag", j + 1, "B", j
  end
  for j in 1:N
    ampo += U/2, "N(N-1)", j
    ampo += -μ, "N", j
  end
  H = MPO(ampo, sites)

  # Initialize the DMRG algorithm with a product state
  psi0 = randomMPS(sites, 1)
  sweeps = Sweeps(max_sweeps)
  setmaxdim!(sweeps, m)
  setcutoff!(sweeps, 1e-12)

  # Perform DMRG sweeps
  energy, psi = dmrg(H, psi0, sweeps)

  # Calculate observables
  println("Ground state energy: ", energy)

  for j in 1:N
    n = real(expect("N", psi, j))
    println("Site $j occupation number: ", n)
  end
end

# Example usage:
N = 10 # Number of lattice sites
t = 1.0 # Hopping amplitude
U = 2.0 # On-site interaction energy
μ = 1.0 # Chemical potential
max_bosons = 5 # Maximum number of bosons per site
max_sweeps = 5 # Number of DMRG sweeps
m = 20 # Number of retained states

bose_hubbard_dmrg(N, t, U, μ, max_bosons, max_sweeps, m)
