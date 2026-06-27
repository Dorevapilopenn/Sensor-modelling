function c_spec = NewtonRaphson(Model, beta, c_tot, c, i)
% Fast Newton-Raphson solver for component/species equilibrium.

persistent cachedModel cachedModelT

if isempty(cachedModel) || ~isequal(cachedModel, Model)
    cachedModel = Model;
    cachedModelT = Model.';
end

c_tot = c_tot(:).';
beta = beta(:).';
c = c(:).';
c_tot(c_tot == 0) = 1e-15;          % numerical difficulties if c_tot=0

for it = 1:100
    % c_spec = beta .* prod(c.^Model, 1), evaluated without repmat/base arrays.
    c_spec = beta .* prod(c.' .^ Model, 1);
    c_tot_calc = c_spec * cachedModelT;
    d = c_tot - c_tot_calc;          % diff actual and calc total conc

    if all(abs(d) < 1e-15)           % return if all diff small
        return
    end

    % J_s(j,k)=sum(Model(j,:).*Model(k,:).*c_spec), vectorized.
    J_s = (Model .* c_spec) * cachedModelT;

    % Solve linear system robustly: use pseudoinverse if Jacobian is ill-conditioned.
    if any(~isfinite(J_s(:)))
        warning('Jacobian contains NaN/Inf at iteration %d (i=%d). Aborting update.', it, i);
        return
    end
    if rcond(J_s) < 1e-12
        delta = d * pinv(J_s);
    else
        delta = d / J_s;
    end

    delta_c = delta .* c;            % equation (2.43), without allocating diag(c)
    c = c + delta_c;

    while any(c <= 0)                % take shift back if conc neg.
        delta_c = 0.5 * delta_c;
        c = c - delta_c;
        if all(abs(delta_c) < 1e-15)
           break
        end
    end
end
