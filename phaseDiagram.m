function phaseDiagram()
% This program is for the generation and the investigation of the Extended 
% Phase Graph (EPG) for a multi-spin-echo pulse sequence. The program is 
% developed based on the paper "Extended phase graphs: Dephasing, RF 
% pulses, and echoes - pure and simple" by M. Weigel 
% (http://doi.org/10.1002/jmri.24619).

% Given the parameters such as flip angle, initial phase and time-position
% of the RF pulses, the complex magnetization is calculated for each of the
% dephasing, rephasing and stored pathways. Each pathway is then drawn on 
% the phase diagram.

% Omer Faruk Oran, May 2016

%% Parameters
% Number of RF pulses including exc RF pulse
RFnum=6;

% Flip angle for the excitation RF pulse
flipAngleExc=pi/2;

% RF phase for the excitation RF pulses.
phiExc=pi/2;

% Flip angle and phase for the refocusing RF pulses.
% turbo-spin-echo
flipAngleRef=ones(1,RFnum)*120/180*pi;
phiRef=zeros(1,RFnum);

% Hyper-echo (use this for testing the program). In hyper-echo, the flip
% angles and the RF phases should be symmetric*
% * Hennig, J., & Scheffler, K. (2001). Hyperechoes. MRM, 46(1), 6–12.
% flipAngleRef=[80 120 180 -120 -80]/180*pi; RFnum=6; % hyper-echo for RFnum=6
% phiRef=[50 120 0 -120 -50]/180*pi-pi/2; % hyper-echo for for RFnum=6

% Default dephasing slope
m0=0.5;

% Crusher Gradients

% no crusher
% crusher = ones(1,(RFnum-1)*2)*m0;

% Regular (equal) turbo-spin-echo crushers
crusher=ones(1,(RFnum-1)*2)*2;

% Crushers for getting approximately only primary echo**. RFnum should be odd.
% The design of crusher gradients is not complete!
% ** Poon, C. S., & Henkelman, R. M. (1992). "Practical T2 quantitation
% for clinical applications". JMRI, 2(5), 541–553.

% crusher=zeros(1,(RFnum-1)*2);
% crusher(1:4:end)=1:(RFnum-1)/2;
% crusher(2:4:end)=1:(RFnum-1)/2;
% crusher(3:4:end)=-(1:(RFnum-1)/2);
% crusher(4:4:end)=-(1:(RFnum-1)/2);
% crusher=fliplr(crusher);


% Width of the crusher regions
ds=0.1;

% Echo positions in time
ect=[0.5 1:RFnum]; % CPMG condition


%% Variables

% Nline is the # of total lines in the phase diagram. Each line represents
% a portion of an echo pathway (FID pathways are neglected). Also each
% refocusing pulse has crusher regions one before and one after, therefore
% the number of lines are increased by the factor of 3. lineMtx defines
% each line in the phase diagram. Columns of lineMtx represents:
% 1: x coordinate of the start point
% 2: y coordinate of the start point
% 3: x coordinate of the end point
% 4: y coordinate of the end point
% 5: Phasing Status: 1: dephasing, 0:stored, -1: dephasing
% 6: Color index
Nline=3*(0.5*(3^(RFnum)-1))-1;
lineMtx=zeros(Nline,6);

% Complex magnetization for each line
M=zeros(Nline,1);

% Number of regions. The crusher regions are counted separately
Nreg=3*RFnum-1;

% Slope of the lines at each region. In crusher regions, The slope is
% increased relative to the applied crusher gradient
mReg=zeros(1,Nreg)+m0;
ind0=ones(1,Nreg);
ind0(1:3:Nreg)=0;
mReg(ind0==1)=[crusher m0];

%% Initial Values

% Start and end points for the first two line
lineMtx(1,1)=ect(1);
lineMtx(1,2)=0;
lineMtx(1,3)=ect(2)-ds;
lineMtx(1,4)=m0*(ect(2)-ect(1)-ds)+lineMtx(1,2);
lineMtx(2,1)=lineMtx(1,3);
lineMtx(2,2)=lineMtx(1,4);
lineMtx(2,3)=ect(2);
lineMtx(2,4)=mReg(2)*ds+lineMtx(2,2);

% Phase Status for the first two line
lineMtx(1:2,5)=1;

% Color index for the first two line
lineMtx(1:2,6)=1;

% Initial magnetization for the first two line
S0=rotT(phiExc,flipAngleExc)*[0; 0; 1];
M(1:2)=S0(1);

%% Calculate pathways and fill the line matrix

u=3; % line index
ci=2; % color index

% the indices of lines which reaches the current refocusing pulse.
reachIndex=2;

for k=1:RFnum-1
    reachIndexAccumulator=[];
    for h=1:length(reachIndex)
        % Fill the line matrix
        mt=mReg(((k-1)*3+1:k*3)+2);
        eny=lineMtx(reachIndex(h),4);
        line0 = lineGenerator(ect(k+1),eny,mt,ect(k+2),ds);
        lineMtx(u:u+8,1:5)=line0;
        lineMtx(u:u+8,6)=ci; % color index
        
        % Calculate the magnetization
        if lineMtx(reachIndex(h),5)==0 %store
            S=[0; 0; M(reachIndex(h))];
        else %dephasing or rephasing
            S=[M(reachIndex(h)); 0; 0];
        end
        S1=rotT(phiRef(k),flipAngleRef(k))*S;
        M(u:3:u+8)=S1(1); % dephasing lines
        M(u+1:3:u+8)=S1(3); % stored lines
        M(u+2:3:u+8)=conj(S1(2)); % rephasing lines
        
        % Increase the index
        u=u+9;
        
        reachIndexAccumulator=[reachIndexAccumulator (u-3:u-1)];
        ci=ci+1;
    end
    reachIndex=reachIndexAccumulator;
end

%% Group related lines and find echo pathways

% Number of pathways
Npath=0.5*(3^RFnum-1); %fids are neglected

% Magnetization for each pathway
pathM=zeros(Npath,1);
pathM(1)=M(1);

% Lines that form the each pathway
pathLine=zeros(3,Npath);
pathLine(1,1)=1;
pathLine(2,1)=2;
pathLine(3,1)=0;

% Group
u=2;
for k=3:9:Nline
    pathM(u)= M(k);
    pathM(u+1)= M(k+1);
    pathM(u+2)= M(k+2);
    
    k0=[k;k+3;k+6];
    pathLine(1:3,u)=k0;
    pathLine(1:3,u+1)=k0+1;
    pathLine(1:3,u+2)=k0+2;
    
    u=u+3;
end

%% Interactive phase diagram

% Default Line Width
lw=2;

hFigure=figure('Position', [300, 300, 1000, 650],'Name', 'Phase Diagram');
movegui(hFigure,'center');

% Colormap
%cmap=colormap(hsv(max(lineMtx(:,6))));
cmap=colormap(lines(max(lineMtx(:,6))));

% String List for the List
formatSpec='%d|||';
formatSpec2='%4.2f';
strList1=num2str((1:length(pathM))',formatSpec);
strList2=num2str(pathM,formatSpec2);
strList=[strList1 strList2];

hList=uicontrol('Style','Listbox','String',strList,...
    'Max',500,'Units','Normalized','Position',[0 0 .11 1]);

%% Crusher plot

% Maximum x value for the plot
xmax=ect(end);

% Axes for the pulse sequence
haxes = axes('Units', 'Normalized','Position', [.15 .75 .80 .17],...
    'XLim', [0 xmax], 'YLim', [-1 1], 'XLimMode', 'manual',...
    'YLimMode', 'manual', 'YTickLabel', '');
title(haxes, 'RF pulses and Crushers (\alpha: Flip angle, \phi:Phase)',...
    'FontSize',12)
hold on

% Draw RF Pulses and Print flip angle and RF phase
fa=[flipAngleExc flipAngleRef]/pi*180;
phi=[phiExc phiRef]/pi*180;
for k=1:RFnum
    line([ect(k); ect(k)],[-1;1],'LineWidth',1,'Color','r');
    text(ect(k)+0.02,0.8,['\alpha=' num2str(fa(k)) '\circ'],'Color','r',...
        'FontSize',12);
    text(ect(k)+0.02,0.5,['\phi=' num2str(phi(k)) '\circ'],'Color','r',...
        'FontSize',12);
end

% Generate and Plot Crusher waveform
N=10000;
crushPlot=zeros(1,N);
c=1;
t=0:xmax/N:xmax-xmax/N;
for k=2:RFnum
    crushPlot(round((ect(k)-ds)/xmax*N):round((ect(k))/xmax*N))=...
        crusher(c)/max(crusher);
    crushPlot(round((ect(k))/xmax*N):round((ect(k)+ds)/xmax*N))=...
        crusher(c+1)/max(crusher);
    c=c+2;
end
plot(t,crushPlot)

%% Interactive Phase Diagram

% Limit y axes
ymax=max(abs(lineMtx(:,4)))*1.1;
% Axes for the phase diagram
haxes2 = axes('Units', 'normalized',...
    'Position', [.15 .050 .80 .62],...
    'XLim', [0 xmax],...
    'YLim', [-ymax ymax]);
title(haxes2, 'Click on a pathway to see its magnetization',...
    'FontSize',14,'Color','r','FontWeight','Bold')

% Draw echo positions and crusher regions
for k=1:RFnum
    line([ect(k); ect(k)],[-ymax;ymax],'LineWidth',1);
    if k>1
        line([ect(k)+ds; ect(k)+ds],[-ymax;ymax],'LineWidth',1,...
            'Color','red','LineStyle','--');
        line([ect(k)-ds; ect(k)-ds],[-ymax;ymax],'LineWidth',1,...
            'Color','red','LineStyle','--');
    end
end
eW=ect(end)-ect(end-1);
line([ect(end-1)+eW/2; ect(end-1)+eW/2],[-ymax;ymax],'LineWidth',1,...
    'Color','magenta','LineStyle','--'); % echo line
line([0; xmax],[0;0],'LineWidth',1); % x-axis

% Draw lines
for k=1:length(lineMtx)
    if lineMtx(k,5)==0
        h=line([lineMtx(k,1); lineMtx(k,3)],[lineMtx(k,2); lineMtx(k,4)],...
            'LineWidth',lw,'Color',cmap(lineMtx(k,6),:),'LineStyle','--');
        set(h,'UserData',k);
        set(h,'buttonDownFcn',@lineClickCallback);
    else
        h=line([lineMtx(k,1); lineMtx(k,3)],[lineMtx(k,2); lineMtx(k,4)],...
            'LineWidth',lw,'Color',cmap(lineMtx(k,6),:));
        set(h,'UserData',k);
        set(h,'buttonDownFcn',@lineClickCallback);
    end
end


%% Internal Functions

%%%% Line Click Callback
    function lineClickCallback(hObject, eventdata)
        
        % get line index
        lineIndex=get(hObject,'UserData');
        % Specify which pathway the line is in
        pathIndex=[find(pathLine(1,:)==lineIndex), find(pathLine(2,:)==lineIndex), ...
            find(pathLine(3,:)==lineIndex)];
        
        % Unhighlight the previous selection
        hline=findobj(gca,'Type','line','-and','LineWidth',lw+1);
        set(hline,'LineWidth',lw);
        
        % Highlight the current selection
        hl1=findobj(gca,'Type','line','-and','UserData',pathLine(1,pathIndex));
        hl2=findobj(gca,'Type','line','-and','UserData',pathLine(2,pathIndex));
        hl3=findobj(gca,'Type','line','-and','UserData',pathLine(3,pathIndex));
        set(hl1,'LineWidth',lw+1);
        set(hl2,'LineWidth',lw+1);
        set(hl3,'LineWidth',lw+1);
        
        % Specify the common lines (pathways) so that they are selected all
        % together in the list. echoIndex gives the indices of the pathways
        % which generate the selected echo
        xs=round(lineMtx(lineIndex,1)*1000);
        ys=round(lineMtx(lineIndex,2)*1000);
        xe=round(lineMtx(lineIndex,3)*1000);
        ye=round(lineMtx(lineIndex,4)*1000);
        f1=round(lineMtx(:,1)*1000)==xs;
        f2=round(lineMtx(:,2)*1000)==ys;
        f3=round(lineMtx(:,3)*1000)==xe;
        f4=round(lineMtx(:,4)*1000)==ye;
        ff=f1+f2+f3+f4;
        sameLineIndex=find(ff>2);
        echoIndex=zeros(length(sameLineIndex),1);
        for k1=1:length(sameLineIndex)
            echoIndex(k1)=[find(pathLine(1,:)==sameLineIndex(k1)), ...
                find(pathLine(2,:)==sameLineIndex(k1)), ...
                find(pathLine(3,:)==sameLineIndex(k1))];
        end
        
        % Select the generating paths in the list
        set(hList,'Value',echoIndex);
        
        % Calculate the total magnetization for the echo
        totalMagnet=sum(pathM(echoIndex));
        hText=findobj(gca,'Type','text');
        delete(hText);
        if imag(totalMagnet)>0
            strText=sprintf(['Total Magnetization in this Echo(Pathway): %4.2f+%4.2fi\n'...
                '(Contributing pathways are highlighted on the left)'],...
                real(totalMagnet),imag(totalMagnet));
        else
            strText=sprintf(['Total Magnetization in this Echo(Pathway): %4.2f%4.2fi\n'...
                '(Contributing pathways are highlighted on the left)'],...
                real(totalMagnet),imag(totalMagnet));
        end
        text(0.1,0.1,strText,'Units','Normalized','FontSize',14,...
            'BackgroundColor','green');
        
    end
end

%% External Functions

%%%% Line Generator Function
function line0 = lineGenerator(stx,sty,m,enx,ds)
% This function obtains the line vectors and end points for one line
% stx: x-coordinate of the start pt
% sty: y-coordinate of the start pt
% m: slope(crusher level) of the 3 states
% enx: x-coordinate of the end point
% ds: spoil width

line0=zeros(9,5);

% 3 lines for rephasing, store, and dephasing for the first region
for k=1:3
    line0(k,1)=stx;
    line0(k,3)=stx+ds;
    if k==1
        line0(k,2)=sty;
        line0(k,4)=m(1)*ds+line0(k,2);
    elseif k==2
        line0(k,2)=sty;
        line0(k,4)=line0(k,2);
    else
        line0(k,2)=-sty;
        line0(k,4)=m(1)*ds+line0(k,2);
    end
    line0(k,5)=(2-k);
end

% 3 lines for rephasing, store, and dephasing for the second region
for k=1:3
    line0(k+3,1)=line0(k,3);
    line0(k+3,2)=line0(k,4);
    line0(k+3,3)=enx-ds;
    if k==1
        line0(k+3,4)=m(2)*(enx-stx-2*ds)+line0(k+3,2);
    elseif k==2
        line0(k+3,4)=line0(k+3,2);
    else
        line0(k+3,4)=m(2)*(enx-stx-2*ds)+line0(k+3,2);
    end
    line0(k+3,5)=(2-k);
end

% 3 lines for rephasing, store, and dephasing for the third region
for k=1:3
    line0(k+6,1)=line0(k+3,3);
    line0(k+6,2)=line0(k+3,4);
    line0(k+6,3)=enx;
    if k==1
        line0(k+6,4)=m(3)*ds+line0(k+6,2);
    elseif k==2
        line0(k+6,4)=line0(k+6,2);
    else
        line0(k+6,4)=m(3)*ds+line0(k+6,2);
    end
    line0(k+6,5)=(2-k);
end
end

%%%% Rotation matrix for the state vector
function f=rotT(phi,alpha)
f=[cos(alpha/2).^2 exp(2i*phi).*sin(alpha/2).^2 -1i*exp(1i*phi).*sin(alpha);
    exp(-2i*phi)*sin(alpha/2).^2 cos(alpha/2).^2 1i*exp(-1i*phi)*sin(alpha);
    -1i/2*exp(-1i*phi)*sin(alpha) 1i/2*exp(1i*phi)*sin(alpha) cos(alpha)];
end