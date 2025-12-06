%% File settings

addpath(genpath("../"))

%% Plot settings

% % Defaults "grid on"
% set(groot,'defaultAxesXGrid','on','defaultAxesXGridMode','manual');
% set(groot,'defaultAxesYGrid','on','defaultAxesYGridMode','manual');
% set(groot,'defaultAxesZGrid','on','defaultAxesZGridMode','manual');
% 
% % Set grid transparency
% set(groot,'defaultAxesGridAlphaMode','manual','defaultAxesGridAlpha',0.25);

% Defaults "grid minor"
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');

% Set grid minor transparency
set(groot,'defaultAxesMinorGridAlphaMode','manual','defaultAxesMinorGridAlpha',0.25);

% IDK
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

% Set font, font style and size
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
set(groot, 'defaultAxesFontSize',14);
set(groot, 'defaultLegendFontSize', 12,'DefaultLegendFontSizeMode','manual');

% Set Figure Position
% set(groot, 'defaultFigurePosition', [470, 360, 700, 430])

% Set Figure Colormap
set(groot, 'defaultFigureColormap', turbo(256));
set(groot, 'defaultSurfaceEdgeAlpha', 0.2);

%Set LineWidth
set(groot, 'defaultLineLineWidth', 1.6);

% Set Figure background colors
set(groot, 'defaultFigureColor', 'w'); % Border color, outside plotted part
set(groot, 'defaultAxesColor', 'w'); % Plotted part backgroud color


