function Hamiltonian = Generate_Hamiltonian(J,U, basis_vectors)
    N=4;
    M=5;
    D = factorial(N+M-1)/(factorial(N)*factorial(M-1));
    
    Hamiltonian = zeros(D,D);
    for u=1:D
        for v=1:D
            u_vect.c = 1;
            u_vect.v = basis_vectors(u,:);
            v_vect.c = 1;
            v_vect.v = basis_vectors(v,:);
           
            element_kinetic = 0;
            for i=1:M
                for j=1:M
                    vect1 = Creation_Operator(j,Annihilation_Operator(i,v_vect));
                    vect2 = Creation_Operator(i,Annihilation_Operator(j,v_vect));
                    element_kinetic = element_kinetic + Bra_Ket(u_vect,vect1) + Bra_Ket(u_vect, vect2);
                end
            end
            
            element_int = 0;
            for i=1:M
                vect1 = Number_Operator(i, Number_Operator(i,v_vect));
                vect2 = Number_Operator(i, v_vect);
                element_int = element_int + Bra_Ket(u_vect, vect1) - Bra_Ket(u_vect,vect2);
            end
            
            element = -J*element_kinetic + U/2*element_int;
            
            Hamiltonian(u,v) = element;
        end
    end
    
end

