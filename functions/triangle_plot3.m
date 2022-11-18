function [] = triangle_plot3()
alph_grid = ... % x1 x2 y1 y2
    [0 1 0 0; 
    1/6 5/6 1/3 1/3; 
    1/3 2/3 2/3 2/3;
    0 1 0 0; 
    1/3 2/3 0 2/3;
    2/3 5/6 0 1/3;
    2/3 1/3 0 2/3;
    1/3 1/6 0 1/3;
    0 1/2 0 1; 
    1/2 1 1 0; ...
    ];

for n = 1:size(alph_grid,1)
    plot3(alph_grid(n, [1 2]), alph_grid(n, [3 4]), .99*[1 1], 'color', [.5 .5 .5], 'linewidth', 1); 
end
text(1/2-.07, 1 + .065, 'AC', 'fontsize', 16)
text(1-.1, -.1, 'CCres',  'fontsize', 16)
text(0-.05, -.1, 'CCunr',  'fontsize', 16)
axis off
xlim([-.1 1.1])
ylim([-.1 1.1])
end
