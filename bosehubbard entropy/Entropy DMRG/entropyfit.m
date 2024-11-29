clear all;

%N = 14; L = 14;
mu = 0; t = 0.2; U=1; V = 0; k = 0; 
dmrg_input.m_warmup = 100;
%dmrg_input.sweep_list = [10,20,25,30,35];


y = [];
w = [];
z = [];
for count=6:20
    L=count;
    N=L;
   conserved_QNum = N;
   model = BoseHubbardChain(L, N, U, mu, t, V, k, BoundCond.periodic, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy1, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;  
   S = real(-trace(M*logm(M))); 
   w(end+1) = S;
   y(end+1) = L;
end
scatter(y,w);