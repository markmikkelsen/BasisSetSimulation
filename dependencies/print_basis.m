function print_basis(BASIS,fullpath_to_basis)
%
xlim_range=[-1 10.0];
%
a=figure;
plot(BASIS.ppm,real(BASIS.specs));legend(BASIS.name,'Location','eastoutside')
set(gca,'xdir','reverse','XGrid','on')
%
text(0.1,0.8,['Echo Time: ',sprintf('%d',BASIS.te)],'Units','normalized')
text(0.1,0.75,[BASIS.seq{1}],'Units','normalized')
text(0.1,0.70,['No. Mets: ',sprintf('%d',BASIS.nMets)],'Units','normalized')
text(0.1,0.65,['LW: ',sprintf('%0.2f',BASIS.linewidth)],'Units','normalized')
text(0.1,0.60,['SpectralW: ',sprintf('%d',BASIS.spectralwidth)],'Units','normalized')
%
ax=gca;
ax.XAxis.MinorTick       = 'on';
ax.XAxis.MinorTickValues = xlim_range(1):0.5:xlim_range(2);
ax.XMinorGrid = 'on';
xlim(xlim_range)
%
exportgraphics(a,fullpath_to_basis, 'Append', false);

for jj=1:BASIS.nMets
    a=figure;
    subplot(2,1,1);
    plot(BASIS.ppm,real(BASIS.specs(:,jj)));
    set(gca,'xdir','reverse','XGrid','on')
    ax=gca;
    ax.XAxis.MinorTick       = 'on';
    ax.XAxis.MinorTickValues = xlim_range(1):0.5:xlim_range(2);
    ax.XMinorGrid = 'on';
    xlim(xlim_range)
    %
    subplot(2,1,2);
    plot(BASIS.ppm,imag(BASIS.specs(:,jj)));
    set(gca,'xdir','reverse','XGrid','on')
    ax=gca;
    ax.XAxis.MinorTick       = 'on';
    ax.XAxis.MinorTickValues = xlim_range(1):0.5:xlim_range(2);
    ax.XMinorGrid = 'on';
    xlim(xlim_range)
    xlim(xlim_range)
    %
    sgtitle(BASIS.name{jj});
    %
    exportgraphics(a,fullpath_to_basis, 'Append', true);
end

end