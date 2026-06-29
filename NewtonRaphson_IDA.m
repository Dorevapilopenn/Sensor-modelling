function c_spec = NewtonRaphson_IDA(beta, c_tot, c, i) %#ok<INUSD>
% Fast safeguarded scalar solver for the IDA model: H, G, D, HD, HG.

beta = beta(:).';
c_tot = c_tot(:).';
c_tot(c_tot <= 0) = 1e-15;

KHD = beta(4);
KHG = beta(5);
Ht = c_tot(1);
Gt = c_tot(2);
Dt = c_tot(3);

lo = 0;
hi = Ht;
H = min(max(c(1), lo), hi);
if H <= 0 || ~isfinite(H)
    H = Ht / max(1 + KHD * Dt + KHG * Gt, 1);
end

tol = max(1e-14, 1e-10 * Ht);
for it = 1:30
    denD = 1 + KHD * H;
    denG = 1 + KHG * H;
    F = H + H * KHD * Dt / denD + H * KHG * Gt / denG - Ht;

    if abs(F) <= tol
        break;
    end

    if F > 0
        hi = H;
    else
        lo = H;
    end

    dF = 1 + KHD * Dt / (denD * denD) + KHG * Gt / (denG * denG);
    Hnew = H - F / dF;
    if Hnew <= lo || Hnew >= hi || ~isfinite(Hnew)
        Hnew = 0.5 * (lo + hi);
    end
    H = Hnew;
end

G = Gt / (1 + KHG * H);
D = Dt / (1 + KHD * H);
HD = KHD * H * D;
HG = KHG * H * G;
c_spec = [H, G, D, HD, HG];
end
