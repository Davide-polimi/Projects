function [] = plot2pdf (inputString)

% plot2pdf è una funzione che permette di salvare un plot in formato pdf
% 
% INPUT: 
% - inputString: è una stringa di testo, semplicente va messo il nome con
% cui si vuole salvare il file
%
% OUTPUT:
% - none
%
% ESEMPIO:
% Se ho un plot e voglio salvarlo con il nome << PLOT_1 >>, mi basterà
% chiamare la funzione plot2pdf e fare:
% plot2pdf ('PLOT_1');
% il file verrà poi salvato (nella cartella in cui è il file di matlab)
% come << PLOT_1.pdf >>
%
% WHY:
% Ho scoperto da poco che overleaf permette di importare le immagini come
% file pdf. Il file pdf ha il vantaggio che se si fa lo zoom non si sgrana
% mai, ideale per avere dei grafici sempre puliti (provare per credere)
%
% AUTORE: Mattia De Marco (+ chatGPT)

title_save = [inputString, '.pdf'];

% Ensure the figure size matches the plot
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
% Export the figure as a PDF
print(gcf, title_save, '-dpdf', '-r0');