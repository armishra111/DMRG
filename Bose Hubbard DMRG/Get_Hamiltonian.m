function Hamiltonian = Get_Hamiltonian(Basis)
    %N is the site count and D the the total number of basis
    [D,N] = size(Basis);
    % the vector is in the form of (constant, vector_elements)
    % index refers to the site number
    function vect = pauli_x(index, basis_vect)
        vect.v = basis_vect.v;
        %rem takes the remainder, so if in state 1 we go to state 0 and
        %vice versa
        vect.v(index) = rem(vect.v(index)+1, 2); 
        vect.c = 1*basis_vect.c;
    end
    
    function vect = pauli_y(index, basis_vect)
        vect.v = basis_vect.v;
        if basis_vect.v(index) == 1
            vect.c = 1i*basis_vect.c;
        else
            vect.c = -1i*basis_vect.c;
        end 
    end

    function vect = pauli_z(index, basis_vect)
        vect.v = basis_vect.v;
        if basis_vect.v(index) == 1
            vect.c = 1*basis_vect.c;
        else
            vect.c = -1*basis_vect.c;
        end 
    end

    function value = Bra_Ket(vect1,vect2)
        if isequal(vect1.v,vect2.v)
            value = vect2.c*vect1.c;
        else
            value = 0;
        end
    end
    j_x = 1;
    j_y = 1;
    j_z = 1;
    h = 0;
    Hamiltonian = zeros(D,D);
    element = 0;
    for i=1: D
        for j=1: D
            flag = 0;
            for k=1: N
                basis_a.v = Basis(i,:);
                basis_a.c = 1;
                basis_b.v = Basis(j,:);
                basis_b.c = 1;
                if k == N
                    k = 1;  
                    flag = 1;
                end
                x_comp = j_x*Bra_Ket(basis_a, pauli_x(k,pauli_x(k+1, basis_b)));
                y_comp = j_y*Bra_Ket(basis_a, pauli_y(k,pauli_y(k+1, basis_b)));
                z_comp = j_z*Bra_Ket(basis_a, pauli_z(k,pauli_z(k+1, basis_b)));
                ex_comp = h*Bra_Ket(basis_a, pauli_z(k, basis_b));
                element = x_comp + y_comp + z_comp + ex_comp;
                Hamiltonian(i,j) = element;
                if flag == 1
                    break;
                end
            end
        end     
    end 
end

