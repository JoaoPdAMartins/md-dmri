function axh = mplot_hist(hist,histpars,axh,opt)

if (nargin < 3), axh = gca; end
if (nargin < 4), opt = []; end

opt = mplot_opt(opt);

if isstruct(hist)
    
    h1 = histogram(axh,hist.b,'BinWidth',histpars.bwidth);
    hold on
    plot([histpars.xmin min(h1.BinEdges)],[0 0],'Color',[0 0 1],'LineWidth',opt.mplot.lw)
    hold on
    plot([max(h1.BinEdges) histpars.xmax],[0 0],'Color',[0 0 1],'LineWidth',opt.mplot.lw)
    hold on
    
    h2 = histogram(axh,hist.g,'BinWidth',histpars.bwidth);
    hold on
    plot([histpars.xmin min(h2.BinEdges)],[0 0],'Color',[0 .7 0],'LineWidth',opt.mplot.lw)
    hold on
    plot([max(h2.BinEdges) histpars.xmax],[0 0],'Color',[0 .7 0],'LineWidth',opt.mplot.lw)
    hold on
    
    h3 = histogram(axh,hist.r,'BinWidth',histpars.bwidth);
    hold on
    plot(axh,[histpars.xmin min(h3.BinEdges)],[0 0],'Color',[1 0 0],'LineWidth',opt.mplot.lw)
    hold on
    plot(axh,[max(h3.BinEdges) histpars.xmax],[0 0],'Color',[1 0 0],'LineWidth',opt.mplot.lw)
    hold off
    
    axis(axh,'tight')
    set([h1 h2 h3],'LineWidth',opt.mplot.lw,'DisplayStyle','stairs',...
        'Normalization',histpars.norm)
    set(h1,'EdgeColor',[0 0 1])
    set(h2,'EdgeColor',[0 .7 0])
    set(h3,'EdgeColor',[1 0 0])
    
    set(axh,'Box','off','TickDir','out')
    set(axh,'FontSize',opt.mplot.fs,'LineWidth',opt.mplot.lw)
else
    
    h1 = histogram(axh,hist);
    hold on

    set(h1,'LineWidth',opt.mplot.lw,'DisplayStyle','stairs',...
        'BinWidth',histpars.bwidth,'Normalization',histpars.norm)
    set(h1,'EdgeColor',[0 0 0])
    
    plot(axh,[histpars.xmin min(h1.BinEdges)],[0 0],'k-','LineWidth',opt.mplot.lw)
    plot(axh,[max(h1.BinEdges) histpars.xmax],[0 0],'k-','LineWidth',opt.mplot.lw)
    hold off
    
    axis(axh,'tight')
    set(axh,'Box','off','TickDir','out')
    set(axh,'FontSize',opt.mplot.fs,'LineWidth',opt.mplot.lw)
end

