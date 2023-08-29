function [sA,sB,sC,sD]=linearizeFEPS(EQ_gnd,M,D,F,PHI,Agrid)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linearize the FEPS (or the NRPS) model
    % x_1 ... x_{ n-1 } are the grounded angles,
    % x_n ... x_{2*n-1} are the frequencies
    % FREISSNER 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nRED=size(Agrid,1);
    %check dimensions
    if(size(D,1)==size(D,2))
        D=diag(D);
        F=diag(F);
        M=diag(M);
    end
    %V: total damping coefficient, ny: part of virtual friction
    V=F+D;
    ny=F./V;

    %calculate Hij=dfi/ddeltaj, Gij = dfi/dthetaj, fi= eq for gen 1
    mn=1/sum((M));
    H=zeros(nRED,nRED);
    G=H;
    for gen=1:1:nRED
        for ogen=1:1:nRED
            if(gen~=ogen)
                h_el=Agrid(gen,ogen)*cos(EQ_gnd(gen)-EQ_gnd(ogen)-PHI(gen,ogen))/M(gen);
                H(gen,ogen)=h_el; %Hij
                H(gen,gen)=H(gen,gen)-h_el;%the sum for the subdiagonal (Hii)
                G(gen,ogen)=mn*M(ogen)/M(gen)*ny(gen)*V(gen); %Gij
            else
                G(gen,gen)=(mn*ny(gen)-1/M(gen))*V(gen); %Gii
            end
        end
    end
    %build system Matrix
    % The matrix is defined as follows
    % x_i = delta_i         for i<=n (note that x_1===0 and will be removed later)
    % x_i = omega_i         for i>n
    % dx_i/dt = dx_i-dx_1   for i<=n
    % dx_(n+i)/dt=Equation  for i>n
    %THE FIRST LINE OF THE MATRIX IS REMOVED LATER
    sA=zeros(2*nRED,2*nRED);
    %fill all elements on the diagonals
    for gen=1:1:nRED
        sA(gen,nRED+gen)=1;
        sA(gen,nRED+1)=-1;
        %d(x_i2)=Equation...
        %put equation coefficients
        sA(nRED+gen,gen)=H(gen,gen);
        sA(nRED+gen,nRED+gen)=G(gen,gen);
        for ogen=1:1:nRED
            if(gen~=ogen)%go through all even rows and add H and G in alternating manner
                sA(nRED+gen,ogen)=H(gen,ogen);
                sA(nRED+gen,nRED+ogen)=G(gen,ogen);
            end
        end
    end
    %Pset is input
    sB=[zeros(nRED,nRED);diag(1./M)];
    
    %y is the vector of all frequencies
    sC=[zeros(nRED,nRED),diag(ones(nRED,1))];
    
    sD=zeros(nRED,nRED);

    %eliminate now the first component of all vectors (delta_1=0)
    sA=sA(2:end,2:end);
    sB=sB(2:end,:);
    sC=sC(:,2:end);
end