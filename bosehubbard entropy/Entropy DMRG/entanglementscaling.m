N=10;L=6; mu = 0.5; t = 0.8; U = 1; V = 0; k = 0; conserved_QNum = NaN;
%model = BoseHubbardChain(L, N, U, mu, t, V, k, BoundCond.periodic, conserved_QNum);
%BoseHubbard_DMRG = DMRG(model);
%mrg_input.m_warmup = 10;
%[Energy, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input);

%M = BoseHubbard_DMRG.reduce_densityM;  

%S = real(-trace(M*logm(M))); 
y = [];
z = []; 
%for ccount=1:20
%   t=ccount/10
for count=4:10
   dmrg_input.m_warmup = 5*count;
   model = BoseHubbardChain(L, N, U, mu, t, V, k, BoundCond.open, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy1, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;  
   S = real(-trace(M*logm(M)));
   z(end+1) = S;
   y(end+1) = dmrg_input.m_warmup;
end
plot(log(y),z);
%save('slist',"y","z")
%end 