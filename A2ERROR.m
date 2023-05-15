clc; close all; clear all;

fileID1 = fopen('M3.txt','r');
fileID2 = fopen('M5.txt','r');
fileID3 = fopen('M7.txt','r');

fileID4 = fopen('uh3.txt','r');
fileID5 = fopen('uh5.txt','r');
fileID6 = fopen('uh7.txt','r');

A = (fscanf(fileID1,'%f',[2 Inf]))';
B = (fscanf(fileID2,'%f',[2 Inf]))';
C = (fscanf(fileID3,'%f',[2 Inf]))';

D = fscanf(fileID4,'%f');
E = fscanf(fileID5,'%f');
F = fscanf(fileID6,'%f');

A(:,3)=D;
B(:,3)=E;
C(:,3)=F;

xyA = A(:,1:2);
xyB = B(:,1:2);
xyC = C(:,1:2);

[common_xz, idxA, idxC1] = intersect(xyA, xyC, 'rows');
[common_yz, idxB, idxC2] = intersect(xyB, xyC, 'rows');

mean(abs(A(idxA,3)-C(idxC1,3))./C(idxC1,3)*100)
mean(abs(B(idxA,3)-C(idxC2,3))./C(idxC2,3)*100)

dlmwrite('ind3.txt', idxA, 'delimiter', ' ');
dlmwrite('ind5.txt', idxB, 'delimiter', ' ');
dlmwrite('ind7.txt', idxC1, 'delimiter', ' ');
