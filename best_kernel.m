function best_phi = best_kernel(S, sigma, target)

    %determine if sigma is stable or unstable

    syms z

    deltaS = -diff(S, z);
    deltaS = matlabFunction(deltaS);


    syms x y real
    H11 = diff(real(deltaS(x+1i*y)), x);
    H12 = diff(real(deltaS(x+1i*y)), y);
    H21 = diff(imag(deltaS(x+1i*y)), x);
    H22 = diff(imag(deltaS(x+1i*y)), y);
    H_fn = matlabFunction([H11, H12; H21, H22]);
    H = H_fn(real(sigma), imag(sigma));
    lambda = eig(H);


    %choose phases to be investigated accordingly
    if lambda <= 0
        phis = linspace(0, pi, 36); %add a small shift to make sure sigma is stable throughout
    else
        phis = linspace(pi, 2*pi, 36);

    end

    %create an array to store the quasipotentials
    Vs_one = zeros(1, 36);
    V_max = 5;
    V_stop = 1000;
    for i = 1:36

       Vs_one(i) = Quasipotential_two_ray_duplicate(S, exp(1i*phis(i)), sigma, target, V_max, 100);


        %if V is maxxed out, increase the limit
       while (Vs_one(i) == V_max) && (V_max < V_stop)

           V_max = 1.5*V_max;
           V_max = min(V_max, V_stop);
           Vs_one(i) = Quasipotential_two_ray_duplicate(S, exp(1i*phis(i)), sigma, target, V_max, 100);


       end




       %V_max is increased/decreased gradually to avoid overestimating V(target)
       if i > 1
           if (Vs_one(i) >= Vs_one(i - 1))

               V_max = max(1.5*Vs_one(i), V_max);
               V_max = min(V_max, V_stop); %cap V_max at 1000

           else

               V_max = 0.75*Vs_one(i);

           end

       end

       disp(i)

    end

    %we repeart the process in the reverse direction, for certainty
    phis = flip(phis);

    Vs_two = zeros(1, 36);
    V_max = 5;
    for i = 1:36

       Vs_two(i) = Quasipotential_two_ray_duplicate(S, exp(1i*phis(i)), sigma, target, V_max, 100);


        %if V is maxxed out, increase the limit
       while (Vs_two(i) == V_max) && (V_max < V_stop)

           V_max = 1.5*V_max;
           V_max = min(V_max, V_stop);
           Vs_two(i) = Quasipotential_two_ray_duplicate(S, exp(1i*phis(i)), sigma, target, V_max, 100);


       end




       %V_max is increased/decreased gradually to avoid overestimating V(target)
       if i > 1
           if (Vs_two(i) >= Vs_two(i - 1))

               V_max = max(1.5*Vs_two(i), V_max);
               V_max = min(V_max, V_stop);

           else

               V_max = Vs_two(i);

           end

       end

    end

    Vs_two = flip(Vs_two);
    phis = flip(phis);
    Vs = min(Vs_one, Vs_two);


    f = figure();

    scatter(phis, Vs, 'x');

    xlabel('Kernel phase')
    ylabel('Quasipotential')
    title('Quasipotential vs arg(K)')
    set(gca, 'YScale', 'linear')

    [best_V, indx] = max(Vs_one);
    best_phi = phis(indx);

end
