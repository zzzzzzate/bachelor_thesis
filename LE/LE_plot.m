% 加载dat数据文件
data = load('lyapunov_results.dat');
% 提取列向量，语法和numpy完全一致
rho = data(:, 1);
psi = data(:, 2);
lyap = data(:, 3);

% 设置画布尺寸
% Position=[100,100,1000,600] 代表 [左偏移,上偏移,宽度,高度]
figure('Position', [100, 100, 1000, 600]);
% 绘制曲线
plot(psi, lyap, 'k.-');

% 设置坐标轴标签
xlabel('\psi');
ylabel('\lambda (Lyapunov exponent)');
% 设置标题
title('Lyapunov Exponent vs Radius');
% 显示网格
grid on;

% 保存图片，指定分辨率150dpi
% print('lyapunov_vs_rho.png', '-dpng', '-r150');
% % 显示图像
% drawnow;