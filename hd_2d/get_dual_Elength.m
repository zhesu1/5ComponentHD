
function Elength = get_dual_Elength(f0, f1)


    if f0<=0 && f1<=0

        Elength = 1;
        
    elseif f0<=0 && f1>0

        Elength = -f0/(f1 - f0);

    elseif f0>0 && f1<=0
        
        Elength = -f1/(f0 - f1);

    else

        Elength = 0;
    end

end
