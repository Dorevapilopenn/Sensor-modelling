function c_spec = NewtonRaphson_m(beta, c_tot, c, i)
% Fast Newton-Raphson solver for the M model: H, G, D, HD, HG, GD.

c_tot = c_tot(:).';
beta = beta(:).';
c = c(:).';
c_tot(c_tot == 0) = 1e-15;          % numerical difficulties if c_tot=0

for it = 1:100
    H = c(1);
    G = c(2);
    D = c(3);

    c_spec = [ ...
        beta(1) * H, ...
        beta(2) * G, ...
        beta(3) * D, ...
        beta(4) * H * D, ...
        beta(5) * H * G, ...
        beta(6) * G * D];

    HD = c_spec(4);
    HG = c_spec(5);
    GD = c_spec(6);

    c_tot_calc = [c_spec(1) + HD + HG, ...
                  c_spec(2) + HG + GD, ...
                  c_spec(3) + HD + GD];
    d = c_tot - c_tot_calc;

    if all(abs(d) < 1e-15)
        return
    end

    J_s = [c_spec(1) + HD + HG, HG,                 HD; ...
           HG,                 c_spec(2) + HG + GD, GD; ...
           HD,                 GD,                 c_spec(3) + HD + GD];

    if any(~isfinite(J_s(:)))
        warning('Jacobian contains NaN/Inf at iteration %d (i=%d). Aborting update.', it, i);
        return
    end
    if rcond(J_s) < 1e-12
        delta = d * pinv(J_s);
    else
        delta = d / J_s;
    end

    delta_c = delta .* c;
    c = c + delta_c;

    while any(c <= 0)
        delta_c = 0.5 * delta_c;
        c = c - delta_c;
        if all(abs(delta_c) < 1e-15)
           break
        end
    end
end

fprintf(1, 'no conv. at C_spec(%i,:)\n', i);
