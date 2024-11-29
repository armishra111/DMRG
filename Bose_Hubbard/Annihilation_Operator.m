function result = Annihilation_Operator(index,vect)
    if vect.v(index) == 0
        vect.c = 0;
    else 
        vect.v(index) = vect.v(index)-1;
        vect.c = vect.c*sqrt(vect.v(index)+1);
    end
    result = vect;
end

