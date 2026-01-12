% 加载dat数据文件
data = load('lyapunov_results.dat');
% 提取列向量
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

%% --- 新增部分：绘制收敛性历史 ---
figure('Position', [150, 150, 1000, 600]);

% 尝试读取收敛数据
if exist('lyapunov_convergence.dat', 'file')
    conv_data = load('lyapunov_convergence.dat');
    % 矩阵尺寸：行=时间步, 列=不同的轨道(不同半径)
    
    % 使用半透明线绘制所有轨道，以便观察整体趋势
    plot(conv_data, 'LineWidth', 0.5);
    
    xlabel('Renormalization Step Index');
    ylabel('\lambda (Lyapunov Exponent)');
    title('Convergence of Lyapunov Exponents over Iterations');
    
    grid on;
    
    % 在图上标注说明
    text(0.05, 0.95, ['Total Orbits: ', num2str(size(conv_data, 2))], ...
        'Units', 'normalized', 'BackgroundColor', 'w');
else
    text(0.5, 0.5, 'lyapunov_convergence.dat not found', ...
        'HorizontalAlignment', 'center');
end

% 保存收敛图
% print('lyapunov_convergence.png', '-dpng', '-r150');