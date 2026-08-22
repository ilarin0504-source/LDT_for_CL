function [x_mesh, y_mesh, V_mesh] = Quasipotential(S, K, sigmas)
%this function uses the ray tracing method to find the Quasipotential for
%the Complex Langevin SDE
    %Evaluate the Freidlin-Wentzel matrix a and vector b
    noise = [sqrt(2)*real(sqrt(K)); sqrt(2)*imag(sqrt(K))];
    a = noise*noise';
    syms z
    deltaS = -K*diff(S, z);
    deltaS = matlabFunction(deltaS);
    b = @(t, w) to_vect(w, deltaS);

    %we will also need the Jacobian
    syms x y real
    Jac11 = diff(real(deltaS(x+1i*y)), x);
    Jac12 = diff(real(deltaS(x+1i*y)), y);
    Jac21 = diff(imag(deltaS(x+1i*y)), x);
    Jac22 = diff(imag(deltaS(x+1i*y)), y);
    Jac = matlabFunction([Jac11, Jac12; Jac21, Jac22]);

    %set the ray parameters
    N = 1000; %number of rays
    epsilon = 0.0001; %initial offset
    t_sim = 5; %the 'flight time' of the test particles before the simulation
            %ends

    V_max = 2; %maximal V that the rays can reach

    %create empty cells to store x, y and V
    %we store info for each critical point separately
    
    N_sigm = numel(sigmas);
    x = cell(N_sigm, 1);
    y = cell(N_sigm, 1);
    V = cell(N_sigm, 1);



    i = 0; %track the number of iterations

    %for each critical point

    for sigma = sigmas
        i = i + 1;

        %create empty arreys to store x, y and V for this critical point
        x_sigm = [];
        y_sigm = [];
        V_sigm = [];


        %the test particles will be launched in all directions
        thetas = linspace(0, 2*pi, N);

        for theta = thetas

            %offset the starting point a bit away from the critical points
            x0 = real(sigma) + epsilon*cos(theta);
            y0 = imag(sigma) + epsilon*sin(theta);

            %use H-J equation to calculate their initial momenta
            n_hat = [cos(theta); sin(theta)];
            b0 = b(0, [x0, y0]);

            %avoid division by 0
            div = n_hat'*a*n_hat;
            if abs(div) < 1e-12
                continue;
            end

            p0 = (-2*(n_hat'*b0)/(n_hat'*a*n_hat))*n_hat;

            %the qusipotential is 0 at sigma
            V0 = 0;

            %to find the later position of the test particle we will
            %integrate the Hamilton equations
            
            q0 = [x0; y0; p0; V0]; %set the initial conditions

            %we stop the integration if any of the values exceed a
            %threshold
            opt = odeset('Events', @(t, q) halt(t, q, V_max));
            [t, q] = ode45(@(t, q) H(t, q, a, b, Jac), [0 t_sim], q0, opt);

            % save the rays
            x_sigm = [x_sigm; q(:, 1)];
            y_sigm = [y_sigm; q(:, 2)];
            V_sigm = [V_sigm; q(:, 5)];



        end

        %save the result for the new critical point
        x{i} = x_sigm;
        y{i} = y_sigm;
        V{i} = V_sigm;
    end

    %interpolate the data:

    %create a grid 
    [x_mesh, y_mesh] = meshgrid(linspace(-2, 2, 2000), linspace(-2, 2, 2000));

    %we interpolate the neighbourhood of each critical point separately, to avoid creating unphysical
    %lines 'connecting' them
    V_mesh = V_max*ones(size(x_mesh)); %preset grid V values to V_max

    for i = 1:N_sigm
        %we interpolate (but don't extrapolate) the values of V
        V_int = scatteredInterpolant(x{i}, y{i}, V{i}, 'linear', 'none'); 
        %hence, if V_int is outside of the region covered during ray
        %tracing, we get NaN
        V_s = V_int(x_mesh, y_mesh);
        %replace it with V_max:
        V_s(isnan(V_s)) = V_max;
        % the 'genuine' V near a critical point is the smallest one
        %(i.e.) not V_max
        V_mesh = min(V_mesh, V_s);
    end


    %visualise the data:

    f = figure();

    %we will plot the Quasipotential alongside with the thimbles

%     [stable_thimbles, unstable_thimbles] = Lefschetz_thimble(S);


    %it may also be helpful to visualise the rays
    
    x_rays = cell2mat(x);
    y_rays = cell2mat(y);
    plot(x_rays, y_rays, 'k.');
    hold on

%     for i = 1:3
%         plot(real(stable_thimbles{i, 1}), imag(stable_thimbles{i, 1}), '-')
%     end


    contour(x_mesh, y_mesh, V_mesh, 10, 'r');
    
    hold off

    axis([-2 2 -2 2])

end

function v = to_vect(w, f)
    in = complex(w(:, 1), w(:, 2));
    out = f(in);
    v = [real(out); imag(out)];
end

function dqdt = H(t, q, a, b, J)
    x = q(1:2);
    p = q(3:4);
    b_n = b(0, x');
    J_n = J(x(1), x(2));

    % use Hamilton equations to find dxdt and dpdt
    dxdt = a*p + b_n;
    dpdt = -(J_n'*p);

    %Finally, we find dVdt and put it all together
    dVdt = 0.5*(p'*a*p);

    dqdt = [dxdt; dpdt; dVdt];
end

function [value, isterminal, direction] = halt(t, q, V_max)
    value = [V_max - q(5); 4 - abs(q(1)); 4 - abs(q(2))];
    isterminal = [1; 1; 1];
    direction = [0; 0; 0];
end