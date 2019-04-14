clear;
close all;

%% Veja o help da função rcosfir para compreender os parâmetros a seguir
fim=5;
T=1;
f=10;
delta=T/f;
alfa1 = 0.25;
alfa2 = 0.5;
alfa3 = 0.75;

%------------------------
%Raíz de cosseno levantado
%------------------------
h1_raiz = rcosfir(alfa1, fim, f, T, 'sqrt');
%Figure 1
figure
eixo = -fim*T:delta:fim*T;
plot(eixo, h1_raiz, 'Color', 'black')
title('Raíz do cosseno levantado')
xlabel('tempo (intervalo de símbolo)')
grid

hold on

h2_raiz = rcosfir(alfa2, fim, f, T, 'sqrt');
plot(eixo,h2_raiz, 'Color', 'red')

h3_raiz = rcosfir(alfa3, fim, f, T, 'sqrt');
plot(eixo, h3_raiz, 'Color', 'blue')
legend('\alpha=0.25', '\alpha=0.5', '\alpha=0.75');

hold off

%------------------------
%Cosseno levantado
%------------------------
h1_sem_raiz = rcosfir(alfa1, fim, f, T);
%Figure 2
figure
eixo = -fim*T:delta:fim*T;
plot(eixo, h1_sem_raiz, 'Color', 'black')
title('Cosseno levantado')
xlabel('tempo (intervalo de símbolo)')
grid

hold on

h2_sem_raiz = rcosfir(alfa2, fim, f, T);
plot(eixo, h2_sem_raiz, 'Color', 'red')

h3_sem_raiz = rcosfir(alfa3, fim, f, T);
plot(eixo, h3_sem_raiz, 'Color', 'blue')
legend('\alpha=0.25', '\alpha=0.5', '\alpha=0.75');

hold off

%------------------------
%Sinal transmitido (Raíz do cosseno levantado)
%------------------------
bits = [1 0 0 1 0]
s = 2*bits-1;
s_up = upsample(s, f);
x1_raiz = conv(s_up, h1_raiz);
%Figure 3
figure
eixo = -fim*T:delta:(fim+(length(bits)-1)+(1-1/f))*T;
plot(eixo, x1_raiz, 'black');
title('Sinal transmitido (Raíz do cosseno levantado)')
xlabel('tempo (intervalo de símbolo)')
xlim([-3 8])
grid

hold on 

x2_raiz = conv(s_up, h2_raiz);
plot(eixo, x2_raiz, 'red');

x3_raiz = conv(s_up, h3_raiz);
plot(eixo, x3_raiz, 'blue');
legend('\alpha=0.25', '\alpha=0.5', '\alpha=0.75');

hold off

%------------------------
%Sinal transmitido (Cosseno levantado)
%------------------------
bits = [1 0 0 1 0]
s = 2*bits-1;
s_up = upsample(s, f);
x1_sem_raiz = conv(s_up, h1_sem_raiz);
%Figure 4
figure
eixo = -fim*T:delta:(fim+(length(bits)-1)+(1-1/f))*T;
plot(eixo, x1_sem_raiz, 'black');
title('Sinal transmitido (Cosseno levantado)')
xlabel('tempo (intervalo de símbolo)')
xlim([-3 8])
grid

hold on

x2_sem_raiz = conv(s_up, h2_sem_raiz);
plot(eixo, x2_sem_raiz, 'red');

x3_sem_raiz = conv(s_up, h3_sem_raiz);
plot(eixo, x3_sem_raiz, 'blue');
legend('\alpha=0.25', '\alpha=0.5', '\alpha=0.75');

hold off

%------------------------
%Resposta em Frequencia Usando a Expressão Exata da Transformada de Fourier (pode ser feito via FFT a partir do pulso no tempo também)
%------------------------
conta = 0;
W = 1;
f1 = (1-alfa2)*W;

for ff = -2*(2*W-f1):0.01:-(2*W-f1)-0.01
    conta = conta + 1;
    H(conta) = 0;
    f(conta) = ff;
end 
for ff=-(2*W-f1):0.01:-f1-0.01
    conta=conta+1;
    H(conta)=(1/(4*W))*(1+cos(pi/(2*W*alfa2)*(abs(ff)-W*(1-alfa2))));
    f(conta)=ff;
end    
for ff=-f1:0.01:f1
    conta=conta+1;
    H(conta)=1/(2*W);
    f(conta)=ff;
end
for ff=f1+0.01:0.01:(2*W-f1)
    conta=conta+1;
    H(conta)=(1/(4*W))*(1+cos(pi/(2*W*alfa2)*(abs(ff)-W*(1-alfa2))));
    f(conta)=ff;
end    
for ff=(2*W-f1)+0.01:0.01:2*(2*W-f1)
    conta=conta+1;
    H(conta)=0;
    f(conta)=ff;
end 

figure
plot(f,2*W*H)
axis([-3 3 0 1.1])
grid
legend('\alpha=0.5');
xlabel('\times R/2 (Hz)');
ylabel('H(f)');

