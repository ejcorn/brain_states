function h = COLORBAR(clim)

clim = abs(clim);
h=colorbar; caxis([-clim clim]); h.Ticks = [-clim 0 clim]; h.TickLabels = [-clim 0 clim];