N = 10; L = 10; mu=0; V = 0; k = 0; conserved_QNum = N;
dmrg_input.m_warmup = 30;
len=10

u_list = (1:len)/len;
t_list = (1:len)/(len);
S_list = [];

for i=1: length(t_list)
    t = t_list(i);
    for j=1: length(u_list)
        tic
        U = u_list(j);
        model = BoseHubbardChain(L, N, U, mu, t, V, k, BoundCond.open, conserved_QNum);
        BoseHubbard_DMRG = DMRG(model);
        [Energy, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input);
        M = BoseHubbard_DMRG.reduce_densityM;
        S = real(-trace(M*logm(M))); 
        fprintf("t=%d,mu=%d,S=%d\n", t,mu,S);
        S_list(len+1-j, i) = S;
        toc
    end
end
heatmap(S_list);