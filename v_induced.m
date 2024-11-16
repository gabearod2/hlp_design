function V_p = v_induced(V_inf, omega, r, NB, CofR, Cd, Cl)
% Using the Newton-Raphson iterative scheme/method
    W0 = 1.0;
    h = 0.01;

    for iter = 1:1000
        f1 = f_of_w(V_inf, omega, r, NB, CofR, Cd, Cl, W0-h);
        f2 = f_of_w(V_inf, omega, r, NB, CofR, Cd, Cl, W0+h);
        fprime = 0.5 * (f2 - f1)/h;
        f3 = f_of_w(V_inf, omega, r, NB, CofR, Cd, Cl, W0);

        if fprime == 0
            iter = 1000;
        else
            W1 = W0 - f3/fprime;
            if abs(W1-W0) < 0.0001 
                iter = 1000;
                W0 = W1;
            end
        end
    end

    V_p = W1;
end 