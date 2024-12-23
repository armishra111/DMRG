using ITensors
using LinearAlgebra

function bose_hubbard_dmrg(N::Int, t::Float64, U::Float64, μ::Float64, max_bosons::Int, max_sweeps::Int, m::Int)
  # ... [The previous bose_hubbard_dmrg function code] ...

  # Perform DMRG sweeps
  energy, psi = dmrg(H, psi0, sweeps)

  # Calculate observables
  println("Ground state energy: ", energy)

  for j in 1:N
    n = real(expect("N", psi, j))
    println("Site $j occupation number: ", n)
  end

  # Calculate entanglement entropy for various subsystem sizes
  subsystem_sizes = 1:N÷2
  entropies = Float64[]
  for x in subsystem_sizes
    region = 1:x
    S = entropy(psi, region)
    push!(entropies, S)
  end

  # Estimate the central charge using the Cardy-Calabrese formula
  L = Float64(N)
  logs = [(log(2L/π * sin(π*Float64(x)/L))) for x in subsystem_sizes]
  A = hcat(logs, ones(length(logs)))
  c_estimate, _ = A \ entropies
  c_estimate *= 6

  println("Estimated central charge: ", c_estimate)
end

# ... [The previous example usage code] ...

bose_hubbard_dmrg(N, t, U, μ, max_bosons, max_sweeps, m)
