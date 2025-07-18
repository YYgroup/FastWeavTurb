clc; clear; close all;

A = readmatrix('Sp.dat');
B = readmatrix('Sp_DNS.dat');

r = A(:,1);
Sp = zeros(10,length(r));
Flat = zeros(10,length(r));
r_DNS = B(:,1);
Sp_DNS = zeros(10,length(r_DNS));
Flat_DNS = zeros(10,length(r_DNS));


for p = 2:4
    Sp(p,:) = A(:,p+1);
    Sp_DNS(p,:) = B(:,p+1);
end
Flat = Sp(4,:)./Sp(2,:).^2;
Flat_DNS = Sp_DNS(4,:)./Sp_DNS(2,:).^2;

figure;
loglog(r,Flat, 'ro', 'MarkerSize', 12, 'LineWidth', 2)
hold on;
loglog(r_DNS,Flat_DNS, 'b-', 'MarkerSize', 12, 'LineWidth', 2)

% 提取绘图数据
%data = [xp; slope_even; K41]';
% 输出到.dat文件
%dlmwrite('Sp_zeta.dat', data, 'delimiter', '\t');    