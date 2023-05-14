clc; close all; clear all;

fileID1 = fopen('xx.txt','r');
fileID2 = fopen('yy.txt','r');


xx = fscanf(fileID1,'%f');
yy = fscanf(fileID2,'%f');

M2(:,1)=xx;
M2(:,2)=yy;

frontr=boundary(xx,yy,0.8); %frontera exterior

% Encontrar los bordes que encierran agujeros o contornos interiores
in1=find (xx>=0.2 & xx<=1.2 & yy==0); %lado inferior
in2=find(yy>=0 & yy<=0.2 & xx==0.2); %lateral izquierdo
in3=find(xx>=0.2 & xx<=1.2 & yy==0.2);%lado superior
in4=find(yy>=0 & yy<=0.2 & xx==1.2); %lateral derecho

n1=length(in1); %número de puntos del lado inferior
n2=length(in2); %número de puntos del lado izquierdo
n3=length(in3); %número de puntos del lado superior
n4=length(in4); %número de puntos del lado derecho


a=in1(2);
in1(2)=[];
in1(n1)=a;

b=in3(2);
in3(2)=[];
in3(n1)=b;

in1=flip(in1);
in4=flip(in4);

in1(1)=[];
in2(1)=[];
in3(1)=[];
in4(1)=[];

n1=length(in1); %número de puntos del lado inferior
n2=length(in2); %número de puntos del lado izquierdo
n3=length(in3); %número de puntos del lado superior
n4=length(in4); %número de puntos del lado derecho

N=n1+n2+n3+n4;

%unificación de los lados en un solo vector
frontd=zeros(N,1); 
for i=1:1:n1
    frontd(i)=in1(i);
end

for i=1:1:n2
    frontd(i+n1)=in2(i);
end

for i=1:1:n3
    frontd(i+n1+n2)=in3(i);
end

for i=1:1:n4
    frontd(i+n1+n2+n3)=in4(i);
end

%Guardado en un .txt de las fronteras para luego abrirlas en Octave
dlmwrite('frontd.txt', frontd, 'delimiter', ' ');
dlmwrite('frontr.txt', frontr, 'delimiter', ' ');
dlmwrite('M5.txt', M2, 'delimiter', ' ');

figure(1)
plot(xx,yy,'x')
hold on
plot(xx(frontr),yy(frontr),"-r","LineWidth",2);
plot(xx(frontd),yy(frontd),'-g',"LineWidth",2);
title('Contorno exterior e interior')
legend('mallado','contorno exterior','contorno interior')
hold on
grid on
grid minor

figure(2)
plot(xx,yy,'x')
hold on
plot(xx(frontr),yy(frontr),'o','color','k',"LineWidth",2);
plot(xx(in1),yy(in1),'o','color','y',"LineWidth",1);
plot(xx(in2),yy(in2),'x','color','b',"LineWidth",1);
plot(xx(in3),yy(in3),'o','color','r',"LineWidth",1);
plot(xx(in4),yy(in4),'x','color','c',"LineWidth",1);
title('Comprobación de puntos de contorno interior y exterior')
legend('mallado','contorno exterior','contorno  inferior','contorno lateral izq','contorno superior','contorno derecho')
hold on
grid on
grid minor