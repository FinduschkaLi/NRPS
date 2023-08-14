function YRed=reduceNetwork(nGEN,YLinesArray3Phase,YLoadsArray3Phase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the reduced admittance matrix YRed from a given 
% set of tie-lines and loads. Output filters are considered loads and
% tie-lines in this model.
% Generator numbers are 1 to nGen
% YLinesArray3Phase is a column array, where each element denotes one line
% admittance: [IDOfNode1,IDOfNode2,Admittance]. numberofNode is
% according to your network.
% YLoadsArray3Phase is a column array, where each elements denotes one load
% admittance: [IDofNode,Admittance]
% in 3-phase systems, multiply the one-line admittance by three.
% FREISSNER 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nFull=max([YLinesArray3Phase(:,1); YLinesArray3Phase(:,2);YLoadsArray3Phase(:,1)]);
    nLines=length(YLinesArray3Phase(:,1));
    nLoads=length(YLoadsArray3Phase(:,1));
    Y=zeros(nFull,nFull);
    %%%%%%%%%%%%
    %add lines
    %%%%%%%%%%%%
    for i=1:1:nLines
        j=YLinesArray3Phase(i,1);
        k=YLinesArray3Phase(i,2);
        Y(j,k)=Y(j,k)-YLinesArray3Phase(i,3);
        Y(k,j)=Y(k,j)-YLinesArray3Phase(i,3);
        Y(j,j)=Y(j,j)+YLinesArray3Phase(i,3);
        Y(k,k)=Y(k,k)+YLinesArray3Phase(i,3);
    end
    %%%%%%%%%%%%
    %sanity check
    %%%%%%%%%%%%
    for i=1:1:nFull
        S=abs(sum(Y(i,:)));
        if(S>1e-13)
            warning(sprintf("G.Error with the sum, value is:%.2f",S));
        end
    end
    %%%%%%%%%%%%
    %add loads
    %%%%%%%%%%%%
    for i=1:1:nLoads
        j=YLoadsArray3Phase(i,1);
        Y(j,j)=Y(j,j)+YLoadsArray3Phase(i,2);
    end
    %%%%%%%%%%%%
    %Reduce network
    %%%%%%%%%%%%    
    Ya=Y(1:nGEN,1:nGEN);
    Yb=Y(1:nGEN,nGEN+1:end);
    Yc=Y(nGEN+1:end,1:nGEN);
    Yd=Y(nGEN+1:end,nGEN+1:end);
    YRed=Ya-Yb*inv(Yd)*Yc;


end


