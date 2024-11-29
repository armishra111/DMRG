%clear all;

N = 14; L = 14; mu = 0; t = 1; V = 0; k = 0; conserved_QNum = N;
dmrg_input.m_warmup = 30;
%dmrg_input.sweep_list = [10,20,25,30,35];


y = [];
w = [];
z = [];
for count=0:30
   U=2+0.1*count;
   model = BoseHubbardChain((L/2), (N/2), U, mu, t, V, k, BoundCond.open, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy1, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;  
   S1 = real(-trace(M*logm(M))); 

   model = BoseHubbardChain((L/2)-1, (N/2)-1, U, mu, t, V, k, BoundCond.open, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy2, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;
   S2 = real(-trace(M*logm(M)));

   c1 = (3/log(cos(pi/L)))*(S1-S2);%closebc
   c2=  2*c1;%open bc
   %z(end+1) = c1;
   w(end+1) = c2;
   y(end+1) = U;
end
plot(y,w)
%save("centralchargeval","y","w","z")