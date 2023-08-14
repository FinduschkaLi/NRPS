function dTET=NK(t,TET,D,P,A,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equation of the Nonlinear Kuramoto model,
    % the inertia-less version of the NRPS model
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate angular difference matrix
    DELTATET=repmat(TET,1,n)-repmat(TET',n,1);
    %calculate model
    dTET=D\(P-A.*sin(DELTATET-PHI)*ones(n,1));
end
