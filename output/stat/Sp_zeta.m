clc; clear; close all;

A = readmatrix('Sp.dat');
[row_num, col_num] = size(A);

%37
i_e = row_num - 4
i_s = i_e - 5
%74

%i_e = row_num-1
%i_s = 5

r = A(i_s:i_e,1);
Sp = zeros(10,i_e-i_s+1);
slope = zeros(1,10);
slope_even = zeros(1,5);

for p = 2:10
    Sp(p,:) = A(i_s:i_e,p+1);

    scaling = polyfit(log(r), log(abs(Sp(p,:))), 1);
    slope(p) = scaling(1);
    if(mod(p,2)==0)
        % 生成拟合直线上的点，用于绘图
        x_fit = linspace(min(log(r)), max(log(r)), 100);
        y_fit = polyval(scaling, x_fit);
        
        % 绘制原始数据点和拟合直线
        figure;
        plot(log(r), log(abs(Sp(p,:))), 'ro', 'DisplayName', 'a', 'MarkerSize', 12);
        hold on;
        plot(x_fit, y_fit, 'b-', 'DisplayName', 'b', 'LineWidth', 2);
        hold off;
        xlabel('x', 'FontSize', 16);
        ylabel('y', 'FontSize', 16);
        legend('FontSize', 14);
        grid on;
    end
end
for p = 1:5
    slope_even(p) = slope(2*p);
    r_sp(p) = slope_even(p)./slope_even(1)./p
end
slope_even
figure;
xp = 2:2:10;
plot(xp,slope_even, 'ro', 'MarkerSize', 12, 'LineWidth', 2)
hold on;
K41 = xp/3;
t_t = 2/3+0.25*(xp-2);
SL94 = xp/9 + 2.*(1-(2/3).^(xp/3));
plot(xp,K41, 'k-', 'LineWidth', 2)
plot(xp,t_t, 'go', 'MarkerSize', 12, 'LineWidth', 2)
plot(xp,SL94, 'b-', 'LineWidth', 2)
xlabel('p', 'FontSize', 16);
ylabel('\zeta_p', 'FontSize', 16);
set(gca, 'FontSize', 14);

figure;
xp = 2:2:10;
plot(xp,slope_even./slope_even(1), 'ro', 'MarkerSize', 12, 'LineWidth', 2)
hold on;
K41 = xp/3;
SL94 = xp/9 + 2.*(1-(2/3).^(xp/3));
plot(xp,K41./K41(1), 'k-', 'LineWidth', 2)
plot(xp,SL94./SL94(1), 'b-', 'LineWidth', 2)
xlabel('p', 'FontSize', 16);
ylabel('\zeta_p/\zeta_2', 'FontSize', 16);
set(gca, 'FontSize', 14);

for p = 1:5
    r_SL94(p) = SL94(p)./SL94(1)./p
    r_sp(p) = slope_even(p)./slope_even(1)./p
end
figure;
plot(xp,r_sp, 'ro', 'MarkerSize', 12, 'LineWidth', 2)
hold on;
plot(xp,r_SL94, 'b-', 'MarkerSize', 12, 'LineWidth', 2)

% 提取绘图数据
data = [xp; slope_even; K41]';

% 输出到.dat文件
dlmwrite('Sp_zeta.dat', data, 'delimiter', '\t');    