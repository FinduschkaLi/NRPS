function drawNetwork(nGen,YLinesArray,tend)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function visualizes the graph of a network
    % based on the tie lines and their impedances
    % loads are not represented
    % 2023 F REISSNER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i_=sqrt(-1);
    dt=0.15;
    nNodes=max(max(YLinesArray(:,1:2)));
    nEdges=length(YLinesArray(:,1));
    %initialize random positions and zero velocity
    pos=rand(nNodes,1)*0.1+rand(nNodes,1)*0.1*i_;
    vel=zeros(nNodes,1);
    %connection and coupling strength matrices
    A=zeros(nNodes,nNodes);%A is used to calculate how strong the coupling is between two nodes
    ADJAM=A; %Adjacency matrix
    edges = ones(2,nEdges)*NaN; %positions of edges as pair of complex numbers
    edgesPN=edges; %edges PN holds the number of the two nodes connected by the edge
    %create edges
    for i=1:1:nEdges
        nA=YLinesArray(i,1);
        nB=YLinesArray(i,2);
        w=YLinesArray(i,3);
        A(nA,nB)=-w;
        A(nB,nA)=-w;
        A(nA,nA)=A(nA,nA)+w;
        A(nB,nB)=A(nB,nB)+w;
        edges(:,i)=[pos(nA,1); pos(nB,1)];
        edgesPN(:,i)=[nA;nB];
        ADJAM(nA,nB)=ADJAM(nA,nB)+1;
        ADJAM(nB,nA)=ADJAM(nB,nA)+1;
    end
    LinkStrength=abs(A)/max(max(abs(A-diag(diag(A)))));
    %initialize figure (if already present, delete it)
    id="cn239LFfAJxmcakl931jlASDM";
    h=findall(groot,'Type','figure','Tag',id);
    close(h)
    f=figure("Name","Network","Tag",id);
    % get a handle for the generator nodes
    h_posG=plot(real(pos(1:nGen)),imag(pos(1:nGen)),'o','LineWidth',3,'Color',[0.8 0.4 0]);
    hold all
    % get a handle for the network nodes
    h_pos=plot(real(pos(nGen+1:end)),imag(pos(nGen+1:end)),'o','LineWidth',3,'Color',[0 0.6 0.7]);
    % get a handle for the edges   
    h_edge=plot(real(edges),imag(edges),'Color',[0.6,0.6,0.6]);
    % add text to each node
    for k=1:1:nNodes
        tx(k)=text(real(pos(k,1))+.01,imag(pos(k,1)),num2str(k));
    end
    %animate***********************************************************
    for t=0:dt:tend
        if(~ishandle(f))
            break;
        end
        %calculate the movement of each node
        dirsAll=repmat(pos,1,nNodes)-transpose(repmat(pos,1,nNodes));
        lensAll=abs(dirsAll);lensAllNZ=lensAll;
        lensAllNZ(lensAll==0)=1;
        dirsAll=dirsAll./lensAllNZ;
        force2=(0.5.*dirsAll./(1+lensAllNZ))*ones(nNodes,1);
        dirsEDG=(repmat(pos,1,nNodes)-transpose(repmat(pos,1,nNodes))).*ADJAM.*LinkStrength;
        forceED=(0.4*dirsEDG.*(lensAll-1))*ones(nNodes,1);
        damping=vel*.25;
        vel=vel+(force2-forceED)*dt-damping*dt;
        pos=pos+vel*dt;
        edges=pos(edgesPN);

        %update positions in plot
        set(h_posG,'XData',real(pos(1:nGen)),'YData',imag(pos(1:nGen)));
        set(h_pos,'XData',real(pos(nGen+1:end)),'YData',imag(pos(nGen+1:end)));
        for k = 1:numel(h_edge)
            set(h_edge(k), 'XData', real(edges(:,k)), 'YData', imag(edges(:,k)));
        end
        mmx=abs(min(real(pos))-max(real(pos)));
        for k=1:nNodes
            tx(k).Position=[real(pos(k,1))+mmx*.02,imag(pos(k,1))];
        end
        %keep even ratios
        xlim([min(real(pos))-.1 max(real(pos))+.1]);
        ylim([min(imag(pos))-.1 max(imag(pos))+.1]);
        axis equal
        pause(.00001);

    end
end