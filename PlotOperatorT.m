clear
clc
close all

I = intervalMatrix([1 0.5;-0.2 0.8], [0.1 0.1; 0.05 0.1]);
G(:,:,1) = [0.2 -0.1; 0.1 0.15];
M = matZonotope([1 0.2; 0.3 0.8], G);

BTIM = intervalMatrix(MatIntTimesMatZon(I,M));
IBM = intMult(I,intervalMatrix(M));
e1 = [1; 0];
e2 = [0; 1];

Isamp = myRandomSampling(matZonotope(I),1000);
Msamp = myRandomSampling(M,1000);
N = 100;
IMe1samp = zeros(2, N*N);
IMe2samp = zeros(2, N*N);

k = 0;

for i = 1:N
    for j = 1:N
        k = k + 1;

        A = Isamp(:,:,i) * Msamp(:,:,j);

        IMe1samp(:,k) = A * e1;
        IMe2samp(:,k) = A * e2;
    end
end

BT_IM = intervalMatrix(MatIntTimesMatZon(I,M));

% figure
% subplot(1,2,1)
% hold on
% plot(polytope(IMe1samp))
% plot(MatIntTimesMatZon(I,M)*e1)
% plot(interval(BTIM)*e1)
% plot(interval(IBM)*e1)
% legend( '$\mathcal{I}\mathcal{M}e_1$', ...
%         '$\mathbb{T}_{\mathcal{I}}(\mathcal{M})e_1$',...
%         '$\mathbb{B}(\mathbb{T}_{\mathcal{I}}(\mathcal{M}))e_1$',...
%         '{\mathcal{I}}\star$\mathbb{B}(\mathcal{M})e_1$','Interpreter','latex')
% subplot(1,2,2)
% hold on
% plot(polytope(IMe2samp))
% plot(MatIntTimesMatZon(I,M)*e2)
% plot(interval(BTIM)*e2)
% plot(interval(IBM)*e2)
% legend( '$\mathcal{I}\mathcal{M}e_2$', ...
%         '$\mathbb{T}_{\mathcal{I}}(\mathcal{M})e_2$',...
%         '$\mathbb{B}(\mathbb{T}_{\mathcal{I}}(\mathcal{M}))e_2$',...
%         '{\mathcal{I}}\star$\mathbb{B}(\mathcal{M})e_2$','Interpreter','latex')

figure

t = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Primo subplot
ax1 = nexttile;
hold(ax1,'on')

plot(polytope(IMe1samp),[1 2], 'LineStyle', '-.', 'EdgeColor', rgb('Black'), 'LineWidth', 2);
plot(MatIntTimesMatZon(I,M)*e1,[1 2], 'LineStyle', '--', 'EdgeColor', rgb('LightGreen'), 'LineWidth', 2);
plot(interval(BTIM)*e1, [1 2], 'EdgeColor', rgb('SeaGreen'), 'LineWidth', 2);
plot(interval(IBM)*e1, [1 2], 'EdgeColor', rgb('Coral'), 'LineWidth', 2);
axis([0.7 1.6 -0.2 0.3])
box on
set(gca,'FontSize',10)

% Secondo subplot
ax2 = nexttile;
hold(ax2,'on')
plot(polytope(IMe2samp),[1 2], 'LineStyle', '-.', 'EdgeColor', rgb('Black'), 'LineWidth', 2);
plot(MatIntTimesMatZon(I,M)*e2,[1 2], 'LineStyle', '--', 'EdgeColor', rgb('LightGreen'), 'LineWidth', 2);
plot(interval(BTIM)*e2, [1 2], 'EdgeColor', rgb('SeaGreen'), 'LineWidth', 2);
plot(interval(IBM)*e2, [1 2], 'EdgeColor', rgb('Coral'), 'LineWidth', 2);
axis([0.2 1 0.3 0.9])
box on
set(gca,'FontSize',10)
set(gca, 'Position',[0.5400    0.1993    0.4050    0.7590])
% Legenda condivisa
% labels = { ...
%     '$\mathcal{I}\mathcal{M}e_i$', ...
%     '$\mathrm{T}_{\mathcal{I}}(\mathcal{M})e_i$', ...
%     '$\mathrm{B}(\mathrm{T}_{\mathcal{I}}(\mathcal{M}))e_i$', ...
%     '$\mathcal{I}\star\mathrm{B}(\mathcal{M})e_i$'};
% lgd = legend(ax1, labels, 'Interpreter','latex','FontSize',15);
% lgd.Layout.Tile = 'south';
% lgd.NumColumns = 2;

lgd = legend(ax1, ...
    {'LEGENDLEGEND1','LEGENDLEGEND2','LEGENDLEGEND3','LEGENDLEGEND4'}, ...
    'Interpreter','none', ...
    'FontSize',10, ...
    'NumColumns',2);

lgd.Layout.Tile = 'south';

