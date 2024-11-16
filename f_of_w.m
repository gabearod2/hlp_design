function f = f_of_w(V_inf, omega, r, NB, CofR, Cd, Cl, w)
    f = 25.1327412287183*w*r/(NB*CofR)-sqrt(1+(omega*r/(V_inf+w))^2)*...
        (Cl*omega*r-Cd*(V_inf + w));
end