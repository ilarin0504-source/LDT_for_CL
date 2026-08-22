function [Y_vals, B, B_err] = boundary_terms(S, K, N_Y)
%this function can be used to calculate the boundary terms of a list of
%observables, given the action S, the (constant) kernel K, and the
%dimension d

    z_symb = sym('z');
    S_symb = S(z_symb);
    deltaS = -K.*gradient(S_symb, z_symb);
    deltaS = matlabFunction(deltaS, 'Vars', {z_symb});

    f = @(t, w) to_vect(w, deltaS);

    %define the observables
    O_sym = [z_symb^6, z_symb^20];
    n_obs = length(O_sym);

    %find the boundary term operator for each observable
    BOs = cell(1, n_obs);
    for i = 1:n_obs
        BOs{i} = inv_langevin_operator(S_symb, O_sym(i), z_symb, K);
    end

    N = 1000;

    w_0 = 4*(rand(N, 2) - 0.5);
    noise = zeros(1, 2);
    noise(1:2:end) = sqrt(2)*real(sqrt(K));
    noise(2:2:end) = sqrt(2)*imag(sqrt(K));

    [dts, w_out] = sde_solver(f, @(t, w) noise, [0 10], 1e-3, w_0, N, 5, 100);
    n_steps = size(w_out, 3);

    %determine Y range from the distribution
    range = zeros(N*n_steps, 1);
    
    for i = 1:n_steps
        w_i = w_out(:, :, i);
        z_i = w_i(:, 1)+1i*w_i(:, 2);
        range((i-1)*N+1 : i*N) = abs(z_i);
    end

    Y_max = min(1.5*max(range), 10);
    Y_vals = linspace(0, Y_max, N_Y);

    %because of the adaptive step size employed by sde_solver_multidim.m,
    %the distribution is biased towards the edge (i.e. the action at the margins of the distribution is larger => 
%   dt is smaller => points are disproportionally likely to be sampled from their)
%in order to remove this bias, the sampled points are weighter by dt
    B = zeros(N_Y, n_obs);
    B_err = zeros(N_Y, n_obs);
    t_tot = sum(dts);
    weights = dts(:)/t_tot;

    B_Y = zeros(n_steps, n_obs);

    for k = 1:N_Y

        for i = 1:n_steps
            
            z_c = w_out(:, 1, i) +1i*w_out(:, 2, i);
            %ignore the points that lie outside the Y range
            outside = abs(z_c) > Y_vals(k);
            for j = 1:n_obs

                BO = BOs{j}(z_c.');
                BO(outside) = 0;
                B_Y(i, j) = mean(BO);

            end
        end
        for j = 1:n_obs
            B(k, j) = dot(weights, B_Y(:, j));
            variance = dot(weights, abs(B_Y(:, j) - B(k, j)).^2);
            B_err(k, j) = 2*sqrt(variance*sum(weights.^2)); 
            %we take the error to be 2 sigma

        end
    end

    %plot the boundary terms for each observable
    for j = 1:n_obs

        f = figure();

        errorbar(Y_vals, abs(B(:, j)), B_err(:, j));
    end

end

function v = to_vect(w, f)

    in = complex(w(:, 1:2:end), w(:, 2:2:end));
    out = f(in.');
    out = out.';
    v = zeros(size(w));
    v(:, 1:2:end) = real(out);
    v(:, 2:2:end) = imag(out);

end

function L_inv = inv_langevin_operator(S_sym, O_sym, z_sym, K_vec)
    % Find the result of the 'inverse langevin' operator acting on a given
    % observable

    gradS = gradient(S_sym, z_sym);
    gradO = gradient(O_sym, z_sym);
    %find the hessian to get second derivatives of O
    hessO = hessian(O_sym, z_sym);

    L_inv_symb = sum(K_vec.'.*(diag(hessO).'-gradS.'.*gradO.'));

    L_inv = matlabFunction(L_inv_symb, 'Vars', {z_sym});

end
