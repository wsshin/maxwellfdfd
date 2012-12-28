function c = b2r(m)
% Color map from blue to white to red.

if nargin < 1, m = size(get(gcf,'colormap'),1); end

% From [0 0 1] to [1 1 1] to [1 0 0];
if mod(m,2) == 0  % color bar depth is even
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else  % color bar depth is odd
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; 1; ones(m1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 
