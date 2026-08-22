function [ts, dts, x] = sde_solver_multidim(f, sigma, t_span, step_size, x_0, N, cutoff, freq, d, checkpoint_file)
    %this function solves an sde equation dx = fdt+sigmadW
    %it implements an adaptive step size to prevent runaway solutions
    %the Euler-Heun method is used
    %optional: pass checkpoint_file to enable periodic saving and resuming
    epsilon = 0.01;
    x_size = 2000;
    block_size = 1000;
    checkpoint_freq = 36000; %save checkpoint every 1 hour

    %resume from checkpoint if one exists
    if nargin >= 10 && ~isempty(checkpoint_file) && isfile(checkpoint_file)
        fprintf('Resuming from checkpoint: %s\n', checkpoint_file);
        cp = load(checkpoint_file);
        x = cp.x;
        dts = cp.dts;
        ts = cp.ts;
        x_n = cp.x_n;
        t_n = cp.t_n;
        counter = cp.counter;
        x_size = cp.x_size;
    else
        x = zeros(N, 2*d, x_size);
        dts = zeros(x_size, 1);
        ts = zeros(x_size, 1);
        x_n = x_0;
        t_n = t_span(1);
        counter = 1;
    end

    last_checkpoint_time = tic;
    while t_n < t_span(2)
        %choose a minimum step size
        drift = f(t_n, x_n);
        norms = sqrt(sum(drift.^2, 2));
        dt_vals = epsilon./norms;
        dt = min([dt_vals; step_size]);
        x_tilda = x_n + dt*f(t_n, x_n);
        v = normrnd(0, 1, [N, d]);
        v = v(:, repelem(1:d, 2));
        x_new = x_n + 0.5*dt*(f(t_n, x_n) + f(t_n + dt, x_tilda)) + ... 
            0.5*(sigma(t_n, x_n) + sigma(t_n + dt, x_tilda))*sqrt(dt).*v;
        t_n = t_n + dt;
        x_n = x_new;
        x_n(:, 2:2:end) = x_n(:, 2:2:end) .* (abs(x_n(:, 2:2:end)) > 1e-300); %if y-values become to small, round them to 0
        % we only record values post-thermalisation (i.e. after the cutoff time)
        if t_n > cutoff
            %we only record every freq^th step
            if mod(counter, freq) == 0
                %if we ran out of the pre-allocated memory, add a new block
                if (counter/freq > x_size)
                    x(:, :, end + block_size) = 0;
                    dts(end + block_size) = 0;
                    ts(end + block_size) = 0;
                    x_size = x_size + block_size;
                end
                x(:, :, round(counter/freq)) = x_n;
                dts(round(counter/freq)) = dt;
                ts(round(counter/freq)) = t_n;
            end
            counter = counter + 1;
        end

        %save checkpoint periodically
        if nargin >= 10 && ~isempty(checkpoint_file) && toc(last_checkpoint_time) > checkpoint_freq
            fprintf('Saving checkpoint at t = %.4f\n', t_n);
            save(checkpoint_file, 'x', 'dts', 'ts', 'x_n', 't_n', 'counter', 'x_size');
            last_checkpoint_time = tic;
        end
    end
    %discard the unused space
    x = x(:, :, 1:floor(counter/freq));
    dts = dts(1:floor(counter/freq));
    ts = ts(1:floor(counter/freq));
end
