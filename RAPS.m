function dXdt=RAPS(t,X,PSI,M,K,D,R,F,Pset,Qset,Vset,Yv_a,Yev,omega_n,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equation of the Reactive Active Power System model,
    % a third order extension of the NRPS model modeling voltage dynamics
    % set F=zeros(n,n) of no virtual friction is to be used
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %retrieve state variables
    TET=X(1:n);
    OMG=X(n+1:n*2);
    E=X(2*n+1:n*3);
    %obtain phasors, amplitudes, angles of V at the PCC
    V_=Yev*(E.*exp(sqrt(-1)*TET));
    V=abs(V_);
    ZET=angle(V_);
    VMAT=repmat(V,1,n).*repmat(V',n,1);
    DELTAZET=repmat(ZET,1,n)-repmat(ZET',n,1);
    %COI-frequency (for virtual friction only)
    omgC=sum(M*OMG)/sum(diag(M));
    %Coupling between voltages at PCC
    B=VMAT.*Yv_a;
    %update state variables
    dOMG=M\(Pset    +D*(omega_n-OMG) + F*(omgC-OMG)   -B.*sin(DELTAZET-PSI)*ones(n,1));
    dE  =K\(Qset    +R*(Vset-V)                       +B.*cos(DELTAZET-PSI)*ones(n,1));
    dTET=OMG;
    dXdt=[dTET;dOMG;dE];
end
