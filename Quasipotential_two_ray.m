function V_target = Quasipotential_two_ray(S, K, sigma, target, V_max, t_sim)
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



    %turn target into a vector

    target = [real(target), imag(target)];

    %set the ray parameters
    N = 360; %number of rays
    epsilon = 0.0001; %initial offset


    acc = 1e-2; %the accuracy required to 'hit' target


    keep = 30; %the number of 'best' minima that are kept



    %the test particles will be launched in all directions
    thetas = linspace(0, 2*pi, N);
    %create a vector containing distances of closest approach for all
    %thetas
    d_coarse = zeros(1, N);
    V_coarse = zeros(1, N);
    score_coarse = zeros(1, N);
    for i = 1:N

        [d_i, V_i] = dist(thetas(i), sigma, a, b, Jac, epsilon, t_sim, V_max, target);
        d_coarse(i) = d_i;
        V_coarse(i) = V_i;

    end

    %identify the local minima
    prev = circshift(d_coarse, [0, 1]);
    next = circshift(d_coarse, [0, -1]);
    mins = (((d_coarse <= prev) & (d_coarse <= next))) | d_coarse < 0.05;
    indx = find(mins);

    %only keep the first few local minima
    [B,I] = sort(V_coarse(indx), 'ascend');
    indx = indx(I(1:min(keep, numel(indx))));

    %refine the results
    dtheta = 2*pi/N;
%create empty arrays to store the refined results

    Vs = [];
    ds = [];
    ths = [];


    for k = indx
        th_current = thetas(k);
        % refine the launch angle to drive the ray onto the target
        th_opt = fminbnd(@(theta) dist(theta, sigma, a, b, Jac, epsilon, t_sim, V_max, target), th_current - dtheta, th_current + dtheta, optimset('TolX', 1e-9));

        [dmin, Vmin] = dist(th_opt, sigma, a, b, Jac, epsilon, t_sim, V_max, target);

        Vs = [Vs, Vmin];
        ds = [ds, dmin];
        ths =[ths, th_opt];
    end


    %in general, it is possible that several rays succeed in reaching the
    %target; in this case, we pick the ray with the smallest
    %quasipotential:
    cands = ds <= acc;
    if any(cands)
        Vs(~cands) = inf;
        [V_target, winner] = min(Vs);
        d_min = ds(winner);
        theta_min = ths(winner);


    else
        %if no ray reached the target, set V at a target to the cutoff
        %value
        V_target = V_max;
        [d_min, winner] = min(ds);
        theta_min = ths(winner);
    end

    %visualise the data:

    %integrate the ray one last time
%     ray = trace_ray(theta_min, sigma, a, b, Jac, epsilon, t_sim, V_max, target);
%
%     figure();
%     plot(real(sigma), imag(sigma), 'bo', 'MarkerFaceColor', 'b');
%     hold on
%     plot(ray(:,1), ray(:,2), 'k.-');
%     plot(target(1), target(2), 'd', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
%     hold off
%     axis equal;
%     grid on
%     legend('critical points', 'winning ray', 'target', 'Location', 'best');
%     disp([V_target, d_min])


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

function [value, isterminal, direction] = halt(t, q, V_max, target)
    dist_target = sqrt((q(1) - target(1))^2 + (q(2) - target(2))^2);
    value = [V_max - q(5); 4 - abs(q(1)); 4 - abs(q(2)); dist_target - 1e-2];
    isterminal = [1; 1; 1; 1];
    direction = [0; 0; 0; 0];
end

function [d, Vc] = dist(theta, sigma, a, b, Jac, epsilon, t_sim, V_max, target)
%find the closest-approach distance of a given ray
%to avoid overestimating the distance, we interpolate the trajectory

%evaluate the trajectory
    q = trace_ray(theta, sigma, a, b, Jac, epsilon, t_sim, V_max, target);

    if isempty(q)
        d = inf;
        score = inf;
        Vc = inf;
        return;
    end

%extract the trajectory and the quasipotential
    x = q(:, 1:2);
    V = q(:, 5);
    n_steps = size(V);

    d = inf;
    Vc = inf;
    score = inf;

    for i = 1:n_steps-1

        pos = x(i, :);
        vec = x(i+1, :) - x(i, :); %find the vector connecting the current point with the consecutive one

        %parametrise the linear segment as pos+v*sec and find s of maximal
        %approach
        s = ((target - pos)*vec')/dot(vec, vec);


        s = max(0, min(1, s));%s must lie between 0 and 1

        closest = pos + s*vec;
        d_current = norm(target - closest);
        if d_current < d
            %update d of closest approach and the corresponding V
            d = d_current;

            Vc = V(i) + s*(V(i+1) - V(i));
        end
    end

    %score combines Vc and d, penalising rays with too high a
    %quasipotential



end

function q = trace_ray(theta, sigma, a, b, Jac, epsilon, t_sim, V_max, target)

    %offset the starting point a bit away from the critical points
    x0 = real(sigma) + epsilon*cos(theta);
    y0 = imag(sigma) + epsilon*sin(theta);

    %use H-J equation to calculate their initial momenta
    n_hat = [cos(theta); sin(theta)];
    b0 = b(0, [x0, y0]);

    %avoid division by 0
    div = n_hat'*a*n_hat;
    if abs(div) < 1e-12
        q = [];
        return
    end

    p0 = (-2*(n_hat'*b0)/(n_hat'*a*n_hat))*n_hat;

    %the qusipotential is 0 at sigma
    V0 = 0;

    %to find the later position of the test particle we will
    %integrate the Hamilton equations

    q0 = [x0; y0; p0; V0]; %set the initial conditions

    %we stop the integration if any of the values exceed a
    %threshold
    opt = odeset('Events', @(t, q) halt(t, q, V_max, target), 'RelTol', 1e-10, 'AbsTol', 1e-12);
    [t, q] = ode45(@(t, q) H(t, q, a, b, Jac), [0 t_sim], q0, opt);

end
