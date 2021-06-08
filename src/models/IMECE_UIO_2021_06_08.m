%% Setup: Script to generate paper figures: 3x3 Example
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Define the plant and model
Am=[-7 2 4; -2 -1 2;-2 2 -1]; %The model with error
C=[0.5 0 1]; % The output matrix
B=[0;0.7;2]; % The input matrix
[rank(ctrb(Am,B)),rank(obsv(Am,C))] %check obsv and ctrb

Lx=[-5]; % True error in the plant

A=Am+B*Lx*C; %Formulate the true plant

%% Check for stable minimum phase transmission zeros (ASD Property)
P1=B*pinv(C*B)*C; 
P2=eye(3)-P1;
S2=null(C);
W2=S2';
A22=W2*P2*Am*W2';
eig(A22) % If two stable poles, system is ASD



sys_id=ss(Am,B,C,0); %ss model of identified model
sys_real=ss(A,B,C,0); % ss model of true plant


[y_id,t_id,x_id]=impulse(sys_id,6.8); %impulse response of identified sys
[y_real,t_real,x_real]=impulse(sys_real,6.8); %impulse of true plant

figure
plot(t_id,y_id,'-k')
hold on
plot(t_real,y_real,'--k')
grid on
legend('Identified Model','True Plant')
xlabel('Time (s)')
ylabel('$y$ Response [arb.]')
%% Controller Synthesis

Fu=[0 0 0 0 0;0 0 1 0 0;0 -4 0 0 0;0 0 0 0 1;0 0 0 -16 0]; %input generator

gamma_e=1*ones(size(C,1),1) %adaptivity parameter

Theta=ones(size(B,2),size(Fu,1)); %theta for input generator

Bbar=padarray(B,[size(Fu,1),0],0,'post'); % Composite input matrix

Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu]; %composite state matrix
eigAbar=eig(Abar); %check stability of composite matrix
Cbar=padarray(C,[0,size(Fu,1)],0,'post'); %composite output matrix

Q1=1*eye(length(Am)); %lqr Q matrix for state
Q2=46.5*eye(length(Fu)); %lqr Q matrix for input
[Klqr,Slqr,elqr]=lqr(Abar',Cbar',blkdiag(Q1,Q2),95,0); %determine optimal observer gains
Klqr=Klqr.';
Kx=Klqr(1:length(Am),:);  %state estimator gains
Ku=Klqr(length(Am)+1:end,:); %input estimator gains
Ac_bar=Abar-Klqr*Cbar; %composite stabalized system
eig(Ac_bar) 



g_step=3; %step gain
g1=2; %sin(2t) gain
g2=0; %cos(2t) gain
g3=0; %sin(4t) gain
g4=4; %cos(4t) gain
g_missing=0; %missing gain

%% Paper plots: x1 x2 ey u eu
sim_out=sim('UIO_IMECE.slx',25);

figure
ax1=subplot(3,1,1);
plot(sim_out.x.Time,squeeze(sim_out.ex.Data(1,:)),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('$e_{x,1}$ response [arb.]')
ax2=subplot(3,1,2);
plot(sim_out.x.Time,squeeze(sim_out.ex.Data(2,:)),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('$e_{x,2}$ response [arb.]')
ax3=subplot(3,1,3);
plot(sim_out.x.Time,squeeze(sim_out.ex.Data(3,:)),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('$e_{x,3}$ response [arb.]')

figure
ax3=subplot(2,1,1);
plot(sim_out.u.Time,squeeze(sim_out.u.Data),'Color', 'r', 'LineStyle', '-','LineWidth',1.5)
grid on
hold on
plot(sim_out.u.Time,squeeze(sim_out.uhat.Data),'Color', 'k', 'LineStyle', '--','LineWidth',0.8)
xlabel('Time (s)')
ylabel('Input')
legend('$u$','$\hat{u}$')
ax4=subplot(2,1,2);
plot(sim_out.eu.Time,squeeze(sim_out.eu.Data),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('Input Error')
legend('$e_u$')


