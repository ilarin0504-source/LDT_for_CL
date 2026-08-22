function [dts, sample_x, sample_y] = Complex_Langevin_Kernelled(S, K)
%here, we assume that the kernel K is a complex number; in general, the
%simulations become far more complex
    %differentiate S
    syms z
    deltaS = -K*diff(S, z);
    deltaS = matlabFunction(deltaS);

    f = @(t, w) to_vect(w, deltaS);

    %set the number of points to be sampled
    N = 1000;


    %noise amplitiude
    noise = 1;

    w_0 = [0.1*(rand(N, 1)-0.5), 0.1*(rand(N, 1)-0.5)];
    [dts, w_out] = sde_solver(f, @(t, w) [noise*sqrt(2)*real(sqrt(K)), noise*sqrt(2)*imag(sqrt(K))], [0 100], 1e-3, w_0, N, 0, 50);
    sample_x = w_out(:, 1, :);
    sample_y = w_out(:, 2, :);
    sample_x = reshape(sample_x, [numel(sample_x), 1]);
    sample_y = reshape(sample_y, [numel(sample_y), 1]);


    [stable_thimbles, unstable_thimbles] = Lefschetz_thimble(S);

    scatter(sample_x, sample_y, 'x', 'k')

    hold on
    for i = 1:3
        plot(real(stable_thimbles{i, 1}), imag(stable_thimbles{i, 1}), '-')
    end

    for i = 1:3
        plot(real(unstable_thimbles{i, 1}), imag(unstable_thimbles{i, 1}), '--')
    end
    hold off

    axis([-2 2 -2 2])
    title('Complex Langevin distribution when S(z) = (-1 + 4i)z^2 + 2z^4, arg(K) = -arg(1+4i)')
    ylabel('Im(z)')
    xlabel('Re(z)')

end

function v = to_vect(w, f)
    in = complex(w(:, 1), w(:, 2));
    out = f(in);
    v = [real(out), imag(out)];
end
