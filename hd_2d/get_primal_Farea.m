
function Farea = get_primal_Farea(f0, f1, f2, f3)

    if f0<=0 && f1<=0 && f2<=0 && f3<=0
        Farea = 1;

    % one vertex inside
    elseif f0<=0 && f1>0 && f2>0 && f3>0
        Farea = (-f0)/(f1 - f0) * (-f0)/(f2 - f0) / 2;

    elseif f0>0 && f1<=0 && f2>0 && f3>0
        Farea = (-f1)/(f0 - f1) * (-f1)/(f3 - f1) / 2;

    elseif f0>0 && f1>0 && f2<=0 && f3>0
        Farea = (-f2)/(f0 - f2) * (-f2)/(f3 - f2) / 2;

    elseif f0>0 && f1>0 && f2>0 && f3<=0
        Farea = (-f3)/(f1 - f3) * (-f3)/(f2 - f3) / 2;

    % two vertices inside
    elseif f0<=0 && f1<=0 && f2>0 && f3>0
        l0 = (-f0)/(f2 - f0);
        l1 = (-f1)/(f3 - f1);

        Farea = (l0 + l1)/2;

    elseif f0<=0 && f1>0 && f2<=0 && f3>0
        l0 = (-f0)/(f1 - f0);
        l1 = (-f2)/(f3 - f2);

        Farea = (l0 + l1)/2;

    elseif f0>0 && f1>0 && f2<=0 && f3<=0
        l0 = (-f2)/(f0 - f2);
        l1 = (-f3)/(f1 - f3);

        Farea = (l0 + l1)/2;


    elseif f0>0 && f1<=0 && f2>0 && f3<=0
        l0 = (-f1)/(f0 - f1);
        l1 = (-f3)/(f2 - f3);

        Farea = (l0 + l1)/2;

    elseif f0<=0 && f1>0 && f2>0 && f3<=0 
        area0 = ((f2*f2)/((f2-f3)*(f2-f0)) + (f1*f1)/((f1-f3)*(f1-f0)))/2;
        Farea = 1 - area0;

    elseif f0>0 && f1<=0 && f2<=0 && f3>0
        area0 = ((f0*f0)/((f0-f1)*(f0-f2)) + (f3*f3)/((f3-f2)*(f3-f1)))/2;
        Farea = 1 - area0;

    % 3 vertices inside
    elseif f0>0 && f1<=0 && f2<=0 && f3<=0

        area0 = (f0*f0)/(2*(f0-f2)*(f0-f1));
        Farea = 1 - area0;

    elseif f0<=0 && f1>0 && f2<=0 && f3<=0

        area0 = (f1*f1)/(2*(f1-f0)*(f1-f3));
        Farea = 1 - area0;

    elseif f0<=0 && f1<=0 && f2>0 && f3<=0

        area0 = (f2*f2)/(2*(f2-f0)*(f2-f3));
        Farea = 1 - area0;

    elseif f0<=0 && f1<=0 && f2<=0 && f3>0

        area0 = (f3*f3)/(2*(f3-f2)*(f3-f1));
        Farea = 1 - area0;

    % no vertices inside
    else
        Farea = 0;
    end
end
