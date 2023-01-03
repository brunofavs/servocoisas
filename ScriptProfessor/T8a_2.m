%T8a_2.m

close all
clear all

% 2.

load ../data_aula9.mat % t = tempo  em segundos
                    % Pref = \delta_DES
                    % P_track = \delta_REAL
figure(1)
plot(t,Pref,'b-')
hold on
plot(t,P_track,'r-')
ylabel('\delta [mm]')
xlabel('Tempo [s]')
title('(a)')
axis([0 5 -10 130])
legend('\delta_{DES}','\delta_{REAL}')
drawnow

figure(2)
plot(t,Pref,'b-')
hold on
plot(t,P_track,'r-')
ylabel('\delta [mm]')
xlabel('Tempo [s]')
title('(b)')
axis([0.9 1.3 -10 130])
legend('\delta_{DES}','\delta_{REAL}')

% Admitindo, por exemplo, que se tinha obtido p1 = p2 = 100 na parte 1.,
% (o que não é necessáriamente verdade):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.a)
% i. A função de transferência G_omega=omega_REAL/omega_DES = 10000/(s+100)^2
%   (é de 2ª ordem);
% ii. A função de transferência G(s) = delta_REAL/omega_REF = G_omega/s =
%   10000/(s*(s+10)^2) (é de 3ª ardem);
% iii. O erro Er(s) da Fig. 3 é Er=delta_DES-delta_REAL = delta_DES -
%   Er*Kp*G => Er=delta_DES/(1+Kp*G);
% iv. Pelo teorema do valor final, erss(t) = lim(t->0) er(t) = lim(s->0)
%   s*Er(s) = lim(s->0) s*delta_DES/(1+Kp*G);
% v. Considerando uma solicitação em degrau, então delta_DES = 1/s, e
% lim(s->0) s*Er(s) = lim(s->0) 1/(1+Kp*G) = 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.b)
s=tf('s');
G = 10000/(s*(s+100)^2)
figure()
rlocus(G) % apesar de o matlab colocar o controlador C no ramo de
% "feedback", ao contrário do nosso caso em que o controlador C está no
% ramo "feedforward", os pólos são os mesmos em ambas as situações

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.c)
% i.G_delta=delta_REAL/delta_DES = Kp*G/(1+Kp*G) =
%   Kp*10000/(s*(s+100)^2+Kp*10000) = 
%   Kp*10000/(S^3+200*s^2+10000*s+Kp*10000) 
syms Kp EPS; % é necessário instalar a "Symbolic Math Toolbox"
ra = routh([1 200 10000 Kp*10000],EPS)
% i. A fórmula da terceira linha assume valores negativos para Kp>200;
% ii. A fórmula da última linha é sempre positiva;
% iii. Logo o sistema é instável para Kp > 200;

figure()
Kp_vec=1:1:500;
for ind = 1:length(Kp_vec)
    Kp = Kp_vec(ind);
    G_theta = Kp*10000/(s^3+200*s^2+10000*s+Kp*10000);
    plot(Kp,max(real(pole(G_theta))),'+') % O pólo "máximo" passa a ter
    % parte real positiva para Kp>200
    hold on
end
plot(Kp_vec,zeros(size(Kp_vec)),'k')
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.d)
Kp_vec = 1:50;
err=[];
figure()
for ind = 1: length(Kp_vec)
    Kp=Kp_vec(ind);
    G_theta = Kp*10000/(s^3+200*s^2+10000*s+Kp*10000);
    y=lsim(G_theta,Pref,t);
    err = [err sum((P_track-y).^2)];
end
plot(Kp_vec,err,'-o')
axis([22 25 3e+4 5e+4])
figure()
plot(t,Pref)
hold on
plot(t,P_track)
Kp=23;
G_theta = Kp*10000/(s^3+200*s^2+10000*s+Kp*10000);
y=lsim(G_theta,Pref,t);
plot(t,y)

figure()
plot(t,Pref)
hold on
plot(t,P_track)
plot(t,y)
axis([0.9 1.3 -10 130])
% NOTA: as diferenças observadas entre os dados experimentais e os do nosso
% modelo poderão eventualmente ser explicadas pelo facto de eu ter admitido
% os valores p1=p2=100 para omega_REAL/omega_REF, o que possível não é
% verdade.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %2.e)
% i.G_delta=delta_REAL/delta_DES = Kp*(s+10)*G/(1+Kp*(s+10)*G) =
%   Kp*(s+10)*10000/(s*(s+100)^2+Kp*(s+10)*10000) = 
%   Kp*(s+10)*10000/(S^3+200*s^2+(10000+Kp*10000)*s+Kp*100000) 
syms Kp EPS; % é necessário instalar a "Symbolic Math Toolbox"
ra = routh([1 200 (10000+Kp*10000) Kp*100000],EPS)
% i. As fórmulas da terceira e quarta linha assumem sempre valores
%   positivos (para k>0, claro)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %2.f)
% i.G_delta=delta_REAL/delta_DES = Kp*(s+z)*G/(1+Kp*(s+z)*G) =
%   Kp*(s+z)*10000/(s*(s+100)^2+Kp*(s+z)*10000) = 
%   Kp*(s+z)*10000/(S^3+200*s^2+(10000+Kp*10000)*s+Kp*z*10000) 
syms Kp z EPS; % é necessário instalar a "Symbolic Math Toolbox"
ra = routh([1 200 (10000+Kp*10000) Kp*z*10000],EPS)
% i. Aqui a minha interpretação é a seguinte:
%   A fórmula da terceira linha assume valores negativos para (note que no
%   neste caso(-p1)*(-p2) = 10000):
%   ((-p1)*(-p2)-50z)*Kp+(-p1)*(-p2)<0 =>
%   Kp<-p1p2/(p1p2-50z)
%   Mas Kp tem que ser positivo, logo temos que ter (pip2-50z)<0
%   ou seja z>(-p1)(-p2)/50, que é menor (50 vezes) do que (-p1)(-p2). Logo,
%   existem valores de z<(-p1)(-p2) para os quais existem também valores de
%   Kp possíveis em que a terceiora linha assume sinal negativo.
