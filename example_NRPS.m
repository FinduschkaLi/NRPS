%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This model shows how to use the NRPS model in a simple 
% 2 generator grid
% 2023 FREISSNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general definitions
i_=sqrt(-1);        %imaginary number
dt=0.004;           %time step of the model
nRED=2;             %number of generators
%nominal values
omega_n=pi*100;
V_n=230;

%VSM configurations*************************************************
Rg =[1 1];
Lg =[0.02 0.02];
Dp=[2 1];
J=[1 0.5];
Dq=[30 15];
Km=[4 2]*1e3;
%Filter parameters of the inverter (only the "CL" part of the LCL filter)
L_filter=7.5e-4;
R_filter=0.5;
C_filter=6.0000e-06;
R_damp = 5;

%Loads**************************************************************
RL=[25];
LL=[0.2];

%Tie-lines**********************************************************
Lline =[ 0.0100];
Rline=[0.5];

% Power references**************************************************
PsetData = [3000 3000]';
QsetData = [1200 1200]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line impedances***************************************************
Ylines=[];
%Add the virtual stators of the VSMs
Ylines=[Ylines; 1 3 3/(Rg(1)+Lg(1)*i_*omega_n)];
Ylines=[Ylines; 2 4 3/(Rg(2)+Lg(2)*i_*omega_n)];
%Add the filters of the inverters
Ylines=[Ylines; 3 5  3/(R_filter(1)+L_filter(1)*i_*omega_n)];
Ylines=[Ylines; 4 6 3/(R_filter(1)+L_filter(1)*i_*omega_n)];
%Add the tie lines
Ylines=[Ylines;5 6 3/(Rline(1)+Lline(1)*i_*omega_n)];


% Load impedances (impedances connected line-ground*****************
%add the capacitor of the LCL filter of the inverters
YLoads =       [ 3 3/(1/(C_filter*omega_n*i_)+R_damp);
                4 3/(1/(C_filter*omega_n*i_)+R_damp);
                ];
%add the actual loads in the grid
YLoads= [YLoads;5 3/RL(1)+3/(i_*omega_n*LL(1));];


% translate VSM parameters to NRPS parameters***********************
M=diag(J)*omega_n;
D=diag(Dp)*omega_n;
R=diag(Dq)*sqrt(2);
K=diag(Km)/omega_n*sqrt(3);
F=zeros(nRED,nRED);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run both the NRPS and the RAPS model for each set of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial values for voltages, angles and frequencies****************
E0=ones(nRED,1)*V_n;
TET0=rand(nRED,1)*2*pi;%randomly distributed resistances
OMG0=ones(nRED,1)*omega_n;
X0=[TET0; OMG0];
%Init empty arrays**************************************************
%Time vector
tspan=0:dt:15;
%generate reduced Admittance matrix (only used for NRPS model)
YRED=reduceNetwork(nRED,Ylines,YLoads);
[A,PHI,P]=generateNRPSParameters(YRED,E0,PsetData,Dp'*omega_n,omega_n);

%run Simulations*****************************************************
[xNRPS, ~]=     RungeK(@(t,X)NRPS(t,X,PHI,M,D,P,A,nRED),tspan,X0);
%get theta, omega, E, V, POUT, QOUT, grounded angles for each timestep
%unpack values for the NRPS model (second order)
TET=xNRPS(1:nRED,:);
OMG=xNRPS(nRED+1:nRED*2,:);
omgCOI=J*OMG./sum(J);
tetCOI=angle(mean(exp(sqrt(-1)*TET)));
GNDA=mod(TET(2:end,:)-TET(1,:)+pi,2*pi)-pi;
%calculate the output powers****************************************
POUT=ones(nRED,length(tspan));
QOUT=POUT;
for i=1:length(tspan)
    POUT(:,i)=getP(PHI,YRED,E0,TET(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figC=0;
figC=figC+1;
f=figure(figC);
f.Name="Omega";
plot(tspan,OMG,':')
ylabel("\omega [rad/s]")
xlabel("time [s]")
grid on;
xlim([tspan(1) tspan(end)])

figC=figC+1;
f=figure(figC);
f.Name="Active Power";
plot(tspan,POUT',':')
ylabel("P [W]")
xlabel("time [s]")
xlim([tspan(1) tspan(end)])
grid on;


figC=figC+1;
f=figure(figC);
f.Name="Grounded Angles";
plot(tspan,GNDA',':')
ylabel("\delta [rad]")
xlabel("time [s]")
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% format plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LineStyle
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    cfig=figHandles(k);
    h_axes=findobj(cfig,'type','axes');
    splots=length(h_axes);
    for p=1:splots
        h_axe=h_axes(p);
        if(isa(h_axe,'matlab.graphics.axis.Axes'))
            h_lines  = findall(h_axe,'Type','line');
            for i=1:1:length(h_lines)
                set(h_lines(i),'LineWidth',3);
            end
        end
    end
end
%% Font Size
sze=17;
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    cfig=figHandles(k);
    h_lines = findall(cfig,'Type','axes');
    isax= strcmp(get(h_lines,'Tag'),'legend');
    h_lines(isax)=[];
    for i=1:1:length(h_lines)
        set(h_lines(i),'FontSize',sze)
        set(h_lines(i), 'FontName', 'Arial')
    end   
    tx=findall(cfig,'Type','text');
    for i=1:1:length(tx)
        set(tx(i),'FontSize',sze);
        set(tx(i),'Interpreter','Tex');
        set(tx(i), 'FontName', 'Arial');
    end
end
%% Sizing
sizeFig=[1000 600];
margins=[.1 .12 .02 .03];
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    cfig=figHandles(k);
    h_axes=findobj(cfig,'type','axes');
    pos=get(cfig,"Position");
    pos(3)=sizeFig(1);
    pos(4)=sizeFig(2);
    set(cfig,"Position",pos);
    pos=get(h_axes(1),"Position");
    pos(1)=margins(1);
    pos(2)=margins(2);
    pos(3)=1-margins(1)-margins(3);
    pos(4)=(1-margins(2)-margins(4));
    
    set(h_axes(1),"Position",pos);
end
%% Legends
figHandles = findall(0,'Type','figure');
for k=1:1:length(figHandles)
    gcf=figHandles(k);
    set (0, 'currentfigure', gcf); 
    if(~isempty(strfind(gcf.Name,"Grounded")))
        hh = legend("$G_2$");
    else
        hh = legend("$G_1$","$G_2$");
    end
    hh.Interpreter = 'latex';
end