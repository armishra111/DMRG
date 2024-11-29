%Entanglement entropy
clear all;

N = 12; L = 12; mu = 0; U = 1; V = 0; k = 0; conserved_QNum = NaN;
dmrg_input.m_warmup = 30;
dmrg_input.sweep_list = [10,20,25,30,35];


y = [];
z = [];
for count=0:10
   t=count/10;
   model = BoseHubbardChain(L/2, N, U, mu, t, V, k, BoundCond.periodic, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy1, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;  
   S1 = real(-trace(M*logm(M))); 

   model = BoseHubbardChain(L/2-1, N, U, mu, t, V, k, BoundCond.periodic, conserved_QNum);
   BoseHubbard_DMRG = DMRG(model);
   [Energy2, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input); 
   M = BoseHubbard_DMRG.reduce_densityM;

   S2 = real(-trace(M*logm(M)));

   c = 6*L^2/(pi^2)*(S1-S2);
   z(end+1) = c;
   %y(end+1) = S;
end
%% 


%%phase diagram
t_list = (1:10)*0.02;
[mu_lists, t_lists] = phase_bound(t_list, BoundCond.open);
scatter(mu_lists, t_lists,5, 'filled')
%[mu_lists, t_lists] = phase_bound(t_list, BoundCond.open, 1);
%scatter(mu_lists, t_lists,5, 'filled')
view([90,-90])

function [mu_lists, t_lists] = phase_bound(t_list, boundCond)
    t_lists = [];
    mu_lists = [];
    base_N = 5;
    base_L = 5;
    mu = 0;
    U = 1;
    V = 0;
    k = 0;
    dmrg_input.m_warmup = 10;
    dmrg_input.sweep_list = [10,20,30];
    
    for j=0:0
        %L = base_L + j*base_L;
        L = base_L;
        N = base_N + j*base_N;
        for i=1: length(t_list)
           t = t_list(i);
           fprintf("\n calculting with t=%d\n", t);
           
           dmrg_input.L = L;
           dmrg_input.target_QNum = N;
           
           model = BoseHubbardChain(L, N, U, mu,t, V, k,boundCond, N); 
           BoseHubbard_DMRG = DMRG(model);
           energy = BoseHubbard_DMRG.fDMRG(dmrg_input);
           %[energy, gs] = model.ExactGsEnergy();

           model = BoseHubbardChain(L, N+1, U, mu,t, V, k, boundCond, N+1);
           BoseHubbard_DMRG = DMRG(model);
           energy_p = BoseHubbard_DMRG.fDMRG(dmrg_input);
           %[energy_p, gs] = model.ExactGsEnergy();

           model = BoseHubbardChain(L, N-1, U, mu,t, V, k, boundCond, N-1);
           BoseHubbard_DMRG = DMRG(model);
           energy_h = BoseHubbard_DMRG.fDMRG(dmrg_input);
           %[energy_h,gs] = model.ExactGsEnergy();

           mu_p = energy_p-energy;
           mu_h = -energy_h+energy;

           if mu_p <= mu_h
              break; 
           end

           t_lists(end+1) = t;
           mu_lists(end+1) = mu_p;
           t_lists(end+1) = t;
           mu_lists(end+1) = mu_h;
        end
    end
    fileArr = [transpose(t_lists) transpose(mu_lists)];
    writematrix(fileArr, "phase_dagram_data.csv");
end