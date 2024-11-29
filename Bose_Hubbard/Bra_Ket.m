function result = Bra_Ket(bra,ket)
    if ket.c == 0
        result = 0;
    elseif isequal(bra.v, ket.v)
        result = ket.c;
    else
        result = 0;
    end
end

