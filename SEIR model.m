b= 0.00003605%birthrate
d= 0.00003605 %death rate 
c= 0.5079%horizontal transmission factor
e=0.1%exposure rate 
a= 0.01190988%death from disease 
g= 0.14285714%recovery rate 
h=0.00555556 %loss of immunity



So=8850 
Eo= 1000
Io= 140 % (0.5-1.5%) initial infected
Ro= 10
No=10000 %total initial population

% syms C

% R0=(e/(e+b))*(C/(g+b+a));
% c=double(vpasolve(R0==0.99, C))%solving for c value based on choosen R0
% real_R0=(e/(e+b))*(c/(g+b+a))


v=[b; d; c; e; a; g; h];
tvalues= 0:100:1000;
[t,y]=ode45(@(t,y) model(t, y,v), tvalues,[So;Eo;Io;Ro;No]);

figure 
subplot(2,2,1)
plot(t,y(:,1),'LineWidth',2)
title('Susceptible')
xlabel('Time')
ylabel('Population')

subplot(2,2,2)
plot(t,y(:,2),'LineWidth',2)
title('Exposed')
xlabel('Time')
ylabel('Population')

subplot(2,2,3)
plot(t,y(:,3),'LineWidth',2)
title('Infected')
xlabel('Time')
ylabel('Population')

subplot(2,2,4)
plot(t,y(:,4),'LineWidth',2)
title('Recovered')
xlabel('Time')
ylabel('Population')

figure 
plot(t, y(:, 1), t, y(:, 2), t, y(:, 3), t, y(:, 4), 'LineWidth',2);
legend('Susceptible', 'Exposed', 'Infected', 'Recovered');
xlabel('Time');
ylabel('Population');
title('SEIR Model')

figure
plot(t,y(:, 3),'LineWidth',2)
xlabel('Time');
ylabel('Infected');
title('Infected Population Over 3 Years');

function [dydt]= model(t,y,v)
%using odefun()
%t is our time span
%y is a vector of our varriable
%v is a vector of our parameters

S=y(1);
E=y(2);
I=y(3);
R=y(4);
N=y(5);

% v is a vector of parameters in order 

b= v(1);
d= v(2);
c= v(3);
e= v(4);
a= v(5);
g= v(6);
h= v(7);


dSdt=  -c*(S*I/N)+b*N-d*S+h*R; 
dEdt= c*(S*I/N)-e*E-d*E; 
dIdt= e*E-d*I-a*I-g*I; 
dRdt= g*I-d*R -h*R;
dNdt= dSdt+dEdt+dIdt+dRdt;
dydt=[dSdt; dEdt; dIdt; dRdt; dNdt];
end 




