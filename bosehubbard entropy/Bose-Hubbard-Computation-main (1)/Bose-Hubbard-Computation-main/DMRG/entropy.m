clear all;

N = 20; L = 20; mu = 0; U = 1; V = 0; k = 0; conserved_QNum = N;
dmrg_input.m_warmup = 10;
dmrg_input.sweep_list = [10,15,20,25,30];


y = [];
z = [];
for count=0:10
   t=count/20;
   model = BoseHubbardChain(L/2, N, U, mu, t, V, k, BoundCond.open, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy1, gnd_state] = BoseHubbard_DMRG.fDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;
   S1 = real(-trace(Mlogm(M))); 

   model = BoseHubbardChain(L/2-1, N, U, mu, t, V, k, BoundCond.open, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy2, gnd_state] = BoseHubbard_DMRG.fDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;

   S2 = real(-trace(Mlogm(M)));

   c = 6L^2/(pi^2)(S1-S2);
   z(end+1) = c;
   %y(end+1) = S;
end
