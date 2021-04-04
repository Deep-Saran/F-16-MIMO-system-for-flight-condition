clear all;clc;
%Matrices of A and B of the model Paramters
A=[ 0.005142  23.0402  -48.8785  -32.0915  0 0 0 0 0;
   -0.000109 -0.526422  0.997184 -0.004425 0 0 0 0 0;
   -0.000337  2.52708  -0.341902  0.000313 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 0 -0.154099 0.082387 -0.998322 0.05376 0;
    0 0 0 0 -19.2246 -0.893601  0.318845 0 0;
    0 0 0 0  2.29583 -0.000888 -0.278676 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 1 0 0];
A(1,:)=A(1,:)/3.28; A(2,1)=A(2,1)/3.28; A(3,1)= A(3,1)/3.28;
B=[  1.585175  1.585175 -1.049275  -1.049275  0 0;
    -0.033078 -0.033078 -0.0558555 -0.0558555 0 0;
    -2.93107  -2.93107  -0.1058865 -0.1058865 0 0;
     0 0 0 0 0 0  ;   
     0.007199 -0.007199  0.0001785 -0.0001785 0.00735  0.021165;
    -6.7916    6.7916   -8.7234     8.7234    0.414519 3.92325;
    -0.752735  0.752735 -0.1342015  0.1342015 1.51008 -1.96651;
     0 0 0 0 0 0;
     0 0 0 0 0 0];
 B(1,:)=B(1,:)*3.28;
%H matrix: Output of Euler angles
H=[ 0 0 0 1 0 0 0 0 0;  
    0 0 0 0 0 0 0 1 0; 
    0 0 0 0 0 0 0 0 1];
D=0;
%Extended Model with tracking error dynamics
IX=eye(3,3);
IE=eye(3,3);
As=[A zeros(9,6); -IX*H zeros(3,6); -H*A zeros(3,3) -IE];
Bs=[B; zeros(3,6); -H*B];
Hs=[ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
Q=1*eye(15,15);
Q(6,6)=1;Q(8,8)=0.5; 
Q(10,10)=1; 
Q(12,12)=100; Q(14,14)=1;
G=0.01*eye(6,6); %minimum value
[K,S,Eig]=lqr(As,Bs,Q,G)
KF=inv(G)*Bs'*S
%%
sys=ss(As,Bs,Hs,D) % Continuous-time system
Ts=0.1; % Sampling time
sys_d=c2d(sys,Ts,'tustin') % Discrete-time system
G_z=tf(sys_d) % Transfer function G(z)
poles=zpk(G_z) % Zeros, poles and gain
% Discrete-time system: State-space model
Ak=sys_d.a, Bk=sys_d.b, Hk=sys_d.c, Dk=sys_d.d;
[Kfeedback,Kn,Eigenvalues]=dlqr(Ak,Bk,Q,G)
%%
x0=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'; % Initial conditions
Ark=[zeros(9,3);IX;zeros(3,3)]*Ts;%Input
N=2000;
t=1:N;
f=1/750;
r=[-0.1*square(pi*f*t,50);-0.1*square(pi*f*t,50);-0.02*square(pi*f*t,50)];

x(:,1)=(Ak-Bk*Kfeedback)*x0+Ark*r(:,1);
for k=1:N
x(:,k+1)=(Ak-Bk*Kfeedback)*x(:,k)+Ark*r(:,k);
u(:,k)=-Kfeedback*x(:,k);
end

k=0:N; 
figure(1)
plot(k,x(4,:),'b',k,x(8,:),'k',k,x(9,:),'r','LineWidth',2);hold on
plot(r(1,:),'--','LineWidth',1);
plot(r(2,:),'g--','LineWidth',1);
plot(r(3,:),'b--','LineWidth',1);
hold off
legend('{\itx}_4k({\itt})','{\itx}_8k({\itt})','{\itx}_9k({\itt})','{\itr1}({\itt})','{\itr2}({\itt})','{\itr3}({\itt})','FontSize',20); legend boxoff;
%%
figure(2)
plot(k,x(4,:),'b','LineWidth',2);hold on
plot(r(1,:),'--','LineWidth',1); hold off
legend('{\itx}_4k({\itt})','{\itr1}({\itt})','FontSize',20); legend boxoff;
figure(3)
plot(k,x(8,:),'k','LineWidth',2);hold on
plot(r(2,:),'g--','LineWidth',1); hold off
legend('{\itx}_8k({\itt})','{\itr2}({\itt})','FontSize',20); legend boxoff;
figure(4)
plot(k,x(9,:),'r','LineWidth',2);hold on
plot(r(3,:),'b--','LineWidth',1);hold off
legend('{\itx}_9k({\itt})','{\itr3}({\itt})','FontSize',20); legend boxoff;
%%
figure(5)
l=0:N-1;
plot(l,u,'LineWidth',1);
legend('{\itu}(1)','{\itu}(2)','{\itu}(3)','{\itu}(4)','{\itu}(5)','{\itu}(6)','FontSize',20); legend boxoff;
