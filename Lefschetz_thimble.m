function [stable_thimbles, unstable_thimbles] = Lefschetz_thimble(S)
    %find the real and imaginary parts of the action
    u = @(x,y)real(S(x+1i*y));
    v = @(x,y)imag(S(x+1i*y));
    %find its critical points
    syms x y real
    dx = diff(u(x, y), x);
    dy = diff(u(x, y), y);
    coordinates = vpasolve([dx==0, dy==0]);
    x_num = double(coordinates.x);
    y_num = double(coordinates.y);
    sigmas = [x_num, y_num];
    %record the number of critical points
    N = length(x_num);

    %find the Hessian
    H = hessian(S(x+1i*y), [x, y]);
    H = matlabFunction(H, "Vars", {x, y});



    %find the Lefschetz thimbles for each critical point
    stable_thimbles = cell(N, 1);
    unstable_thimbles = cell(N, 1);
    syms z
    deltaS = conj(diff(S, z));
    deltaS = matlabFunction(deltaS);
    for i = 1:N
        %find the eigenvectors
        mat = H(sigmas(i, 1), sigmas(i, 2));
        [V, D] = eig(mat);
        %we are interested in the stable manifold
        [d, ind] = sort(real(diag(D)));
        Vs = V(:, ind);
        %use the hessian to construct an initial point
        epsilon = 0.001;
        z_naught = sigmas(i, 1) + 1i*sigmas(i, 2) + epsilon*(Vs(1, 1) + Vs(2, 1));
        [t_plus, q_plus] = ode45(@(t, z) deltaS(z), [0, 20000], z_naught);
        z_naught = sigmas(i, 1) + 1i*sigmas(i, 2) - epsilon*(Vs(1, 1) + Vs(2, 1));
        [t_minus, q_minus] = ode45(@(t, z) deltaS(z), [0, 20000], z_naught);
        stable_thimbles{i, 1} = [flipud(q_plus); q_minus]; %flip the first array so that it represents a continious path
        %repeat the procedure to find the unstable thimbles
        z_naught = sigmas(i, 1) + 1i*sigmas(i, 2) + epsilon*(Vs(1, 2) + Vs(2, 2));
        [t_plus, q_plus] = ode45(@(t, z) -deltaS(z), [0, 2], z_naught);
        z_naught = sigmas(i, 1) + 1i*sigmas(i, 2) - epsilon*(Vs(1, 2) + Vs(2, 2));
        [t_minus, q_minus] = ode45(@(t, z) -deltaS(z), [0, 2], z_naught);
        unstable_thimbles{i, 1} = [flipud(q_plus); q_minus];
    end


    %Lastly, we visualize the thimble structure

    f = figure();
    hold on
    for i = 1:N
        plot(real(stable_thimbles{i, 1}), imag(stable_thimbles{i, 1}), '-')
        %plot the relevant critical point
        plot(sigmas(i, 1), sigmas(i, 2), 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
        plot(real(unstable_thimbles{i, 1}), imag(unstable_thimbles{i, 1}),'--') 
        xlabel('Re(z)')
        ylabel('Im(z)')
        title('Lefschetz thimbles for S(z) = (1+4i)z^2 + 2z^4 passing through z_3')
    end
    v = -2:0.1:2;
    plot(v, 1.61*v+1.25)
%     %plot two fictious lines for the legend
%     fic_1 = plot(NaN, NaN, '-k');
%     fic_2 = plot(NaN, NaN, '--k');
    hold off
    axis([-2 2 -2 2])
%     %legend([fic_1, fic_2], {'Stable Thimbles', 'Unstable Thimbles'})
%     %plot the contours of Im(S) = const
%     for i = 1:N
%         f1 = figure();
%         fimplicit(@(x, y) v(x, y) - v(sigmas(i, 1), sigmas(i, 2)))
%         xlabel('Re(z)')
%         ylabel('Im(z)')
%         title(['Curves satisfying Im(S(z)) = Im(S(z_))' i+1])
%         axis([-2 2 -2 2])
%     end
end