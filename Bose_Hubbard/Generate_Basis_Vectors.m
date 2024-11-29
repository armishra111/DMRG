function basis_vectors = Generate_Basis_Vectors(N,M)
%Generates the basis vectors given N=number of atoms and M=number of sites
    D = factorial(N+M-1)/(factorial(N)*factorial(M-1));
    basis_vectors = zeros(D,M);
    basis_vectors(1,1) = N;
    current_row = 1;
    while basis_vectors(current_row, M) ~= N
        k = 0;
        for i = M-1: -1: 1
            if basis_vectors(current_row, i) > 0
                k = i;
                break;
            end 
        end 
        
        next_row = current_row + 1;
        for i = 1: M
            if 1 <= i && i <= k-1
                basis_vectors(next_row,i) = basis_vectors(current_row,i);
            elseif i >= k+2
                basis_vectors(next_row,i) = 0;
            end
        end
        basis_vectors(next_row, k) = basis_vectors(current_row, k) - 1;
        basis_vectors(next_row, k+1) = N - sum(basis_vectors(next_row, 1:k));
        current_row = next_row;
    
    end
 
end 
