function result = Creation_Operator(index,vect)
    vect.v(index) = vect.v(index) + 1;
    vect.c = vect.c*sqrt(vect.v(index));
    result = vect;
end

