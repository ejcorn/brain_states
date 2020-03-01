function txt = LABELROUND2(x)

% round real number x to 2 significant digits and convert to string
% used for making plot titles or overlays with stats

txt = num2str(round(x,2,'significant'));
