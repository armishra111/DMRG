N = 10; L = 10; mu = 0.5; t = 0.05; U = 1; V = 0; k = 0; conserved_QNum = NaN;
model = BoseHubbardChain(L, N, U, mu, t, V, k, BoundCond.periodic, conserved_QNum);
BoseHubbard_DMRG = DMRG(model);
dmrg_input.m_warmup = 20;
[Energy, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input);

M = BoseHubbard_DMRG.reduce_densityM;  
S = real(-trace(M*logm(M))); 


