function c_spec=NewtonRaphson(Model, beta, c_tot, c,i)
ncomp=length(c_tot);			% number of components
nspec=length(beta);				% number of species
c_tot(c_tot==0)=1e-15;			% numerical difficulties if c_tot=0

it=0;
while it<=99
    it=it+1;
    % Ensure c is a row vector
    c = reshape(c,1,[]);

    % Compute species concentrations robustly, avoid 0^0 which gives NaN
    base = repmat(c',1,nspec);
    % Wherever exponent (Model) is 0 the contribution should be 1 (not 0^0)
    base(Model==0 & base==0) = 1;
    c_spec    = beta .* prod(base.^Model, 1); % species conc
    c_tot_calc=sum(Model.*repmat(c_spec,ncomp,1),2)'; %comp ctot calc
    d         =c_tot-c_tot_calc;	% diff actual and calc total conc
    
    if all(abs(d) <1e-15)			% return if all diff small
        return
    end
    
    for j=1:ncomp					% Jacobian (J_s=J*)
        for k=j:ncomp
            J_s(j,k)=sum(Model(j,:).*Model(k,:).*c_spec);
            J_s(k,j)=J_s(j,k);      % J_s is symmetric
        end
    end
    % Solve linear system robustly: use pseudoinverse if Jacobian is ill-conditioned
    if any(isnan(J_s(:))) || any(isinf(J_s(:)))
        warning('Jacobian contains NaN/Inf at iteration %d (i=%d). Aborting update.', it, i);
        return
    end
    if rcond(J_s) < 1e-12
        delta = d * pinv(J_s);
    else
        delta = d / J_s;
    end
    delta_c = delta * diag(c);        % equation (2.43)
    c=c+delta_c;
    
    while any(c <= 0)				% take shift back if conc neg.
        delta_c=0.5*delta_c;
        c=c-delta_c;
        if all(abs(delta_c)<1e-15)
           break
        end
	  end
end 

if it>99; fprintf(1,'no conv. at C_spec(%i,:)\n',i); end