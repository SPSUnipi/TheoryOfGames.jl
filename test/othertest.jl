using SankeyPlots

src =     [1,     1,   1,   1,   2,   2, 2,    3,   4,   5, 9, 4]
dst =     [6,     3,   7,   4,   3,   7, 9,    7,   8,   8, 4, 1]
weights = [0.1, 0.3, 0.5, 0.5, 0.2, 2.8, 1, 0.45, 4.5, 3.3, 2., 2.]

sankey(src, dst, weights)