function [dts, x] = sde_solver(f, sigma, t_span, step_size, x_0, N, cutoff, freq)
    %this function solves an sde equation dx = fdt+sigmadW
    %it implements an adaptive step size to prevent runaway solutions
    %the Euler-Heun method is used
    epsilon = 0.01;
    x_size = 200;
    block_size = 100;
    x = zeros(N, 2, x_size);
    dts = zeros(x_size, 1);
    x_n = x_0;
    t_n = t_span(1);
    counter = 1;
    while t_n < t_span(2)
        %choose a minimum step size
        drift = f(t_n, x_n);
        norms = sqrt(sum(drift.^2, 2));
        dt_vals = epsilon./norms;
        dt = min([dt_vals; step_size]);
        x_tilda = x_n + dt*f(t_n, x_n);
        v = normrnd(0, 1, [N, 1]);
        x_new = x_n + 0.5*dt*(f(t_n, x_n) + f(t_n + dt, x_tilda)) + ... 
            0.5*(sigma(t_n, x_n) + sigma(t_n + dt, x_tilda))*sqrt(dt).*[v, v];
        t_n = t_n + dt;
        x_n = x_new;
        % we only record values post-thermalisation (i.e. after the cutoff time)
        if t_n > cutoff
            %we only record every freq^th step
            if mod(counter, freq) == 0
                %if we ran out of the pre-allocated memory, add a new block
                if (counter/freq > x_size)
                    x(:, :, end + block_size) = 0;
                    dts(end + block_size) = 0;
                    x_size = x_size + block_size;
                end
                x(:, :, round(counter/freq)) = x_new;
                dts(round(counter/freq)) = dt;
            end
            counter = counter + 1;
        end
    end
    %discard the unused space
    x = x(:, :, 1:floor(counter/freq));
    dts = dts(1:floor(counter/freq));
end
