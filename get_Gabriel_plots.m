function get_Gabriel_plots(model,par)

% Times at which the data are extracted
tpoint = [0 3.3e-3, 6.7e-3, 10.0e-3, 13.3e-3, 16.7e-3]; 

t1 = mpheval(model,'t','edim',0,'selection',1,'dataset','dset3','solnum','all','dataonly','on');
theta = linspace(0,2*pi,181);

% For images at plane y = 0
rcoords = linspace(0,1.5*par.rcell,81);
zcoords = linspace(-1.5*par.rcell,1.5*par.rcell,161);
[rgrid0,zgrid0] = meshgrid(rcoords,zcoords);
rgrid = [-rgrid0(:,end:-1:1), rgrid0];
zgrid = [zgrid0, zgrid0];
% For convolved image
ycoords = linspace(0,0.5*par.y_opt,61);
[xgrid1,ygrid1,zgrid1] = meshgrid(rcoords,ycoords,zcoords);
factor = exp(-0.5*ycoords.^2/par.y_opt^2);
factor = repmat(factor',1,length(rcoords),length(zcoords));
c12conv_line = [];

figure; hold on
px = 150; % Side length of the subplot in px
sepx = 10; % Separation between columns in px
sepy = 10; % Separation between rows in px

set(gcf,'Position',[80  80   1382  689])
for i = 1:length(tpoint)
    idx = find(t1 >= tpoint(i),1,'first');
    [c10, c20, c120] = mphinterp(model,{'c1','c2','c12'},'coord',[rgrid0(:), zgrid0(:)]','dataset','dset3','solnum',idx);
    c12conv0 = mphinterp(model,'c12','coord',[ygrid1(:),xgrid1(:),zgrid1(:)]','dataset','rev1','solnum',idx);

    c10 = reshape(c10,length(zcoords),length(rcoords));
    c20 = reshape(c20,length(zcoords),length(rcoords));
    c120 = reshape(c120,length(zcoords),length(rcoords));
    c12conv0 = reshape(c12conv0,length(ycoords),length(rcoords),length(zcoords));
    c12conv1 = 2*squeeze(sum(c12conv0.*factor,1).*diff(ycoords(1:2)))';
    c1 = [c10(:,end:-1:1), c10];
    c2 = [c20(:,end:-1:1), c20];
    c12 = [c120(:,end:-1:1), c120];
    c12conv = [c12conv1(:,end:-1:1), c12conv1];
    c12conv_line = [c12conv_line c12conv1(:,1)];
    
    s1 = subplot(4,length(tpoint),i); hold on
    set(gca,'Units','pixels','Position',[(i-1)*(px+sepx)+20 3*(px+sepy)+20 px px])
    surf(zgrid',rgrid',c1')
    plot3(par.rcell*cos(theta),par.rcell*sin(theta),max(max(c1))*ones(size(theta)),'w','LineWidth',2)
    view(2);  colormap jet; shading FLAT
        caxis([0 12]);
    box on; axis image; axis off
    if mod(i,length(tpoint)) == 0
        s1Pos = get(s1,'position');
        hb = colorbar('location','eastoutside');
        set(s1,'position',s1Pos);
    end
    title([num2str(tpoint(i),'%.1e'),' s'])
    
    s2 = subplot(4,length(tpoint),i + length(tpoint)); hold on
    set(gca,'Units','pixels','Position',[(i-1)*(px+sepx)+20 2*(px+sepy)+20 px px])
    surf(zgrid',rgrid',c2')
    plot3(par.rcell*cos(theta),par.rcell*sin(theta),max(max(c1))*ones(size(theta)),'w','LineWidth',2)
    view(2); colormap jet; shading FLAT
    caxis([0 2.3e-3]);
    box on; axis image; axis off
    if mod(i,length(tpoint)) == 0
        s2Pos = get(s2,'position');
        hb = colorbar('location','eastoutside');
        set(s2,'position',s2Pos);
    end
    
    s3 = subplot(4,length(tpoint),i + 2*length(tpoint)); hold on
    set(gca,'Units','pixels','Position',[(i-1)*(px+sepx)+20 1*(px+sepy)+20 px px])
    surf(zgrid',rgrid',c12')
    plot3(par.rcell*cos(theta),par.rcell*sin(theta),max(max(c1))*ones(size(theta)),'w','LineWidth',2)
    view(2); colormap jet; shading FLAT
    caxis([0 2.3e-3]);
    box on; axis image; axis off
    if mod(i,length(tpoint)) == 0
        s3Pos = get(s3,'position');
        hb = colorbar('location','eastoutside');
        set(s3,'position',s3Pos);
    end

    s4 = subplot(4,length(tpoint),i + 3*length(tpoint)); hold on
    set(gca,'Units','pixels','Position',[(i-1)*(px+sepx)+20 0*(px+sepy)+20 px px])
    surf(zgrid',rgrid',c12conv')
    plot3(par.rcell*cos(theta),par.rcell*sin(theta),max(max(c1))*ones(size(theta)),'w','LineWidth',2)
    view(2); colormap jet; shading FLAT; 
    caxis([0 35e-9]);
    box on; axis image; axis off
    if mod(i,length(tpoint)) == 0
        s4Pos = get(s4,'position');
        hb = colorbar('location','eastoutside');
        set(s4,'position',s4Pos);
    end
end

%% 1D concentration profile along the central axis
data = mpheval(model,{'c1','c2','c12'},'edim',1,'selection',[1 3 5],'dataset','dset3','solnum','all','dataonly','off','refine',2);
[z,I] = sort(data.p(2,:)');
c1a1 = data.d1(:,I);
c2a1 = data.d2(:,I);
c12a1 = data.d3(:,I);
c1a = [c1a1];
c2a = [c2a1];
c12a = [c12a1];

idx_t = zeros(size(tpoint));
for i = 1:length(tpoint)
    idx_t(i) = find(t1 >= tpoint(i), 1, 'first');
end

figure; 
set(gcf,'Position',[227         160        1182        1009])
subplot(2,2,1); hold on; box on
plot(z*1e6,c1a(idx_t,:),'LineWidth',1.5)
xlabel('z (\mum)'); ylabel('[c1] (mM)')
xlim([-12 12]); ylim([0 1.1*max(max(c1a))])
% legend('  0.0 ms','  3.3 ms','  6.7 ms','10.0 ms','13.3 ms','16.7 ms')
legend(num2str(tpoint','%.1e'),'Location','NorthEast')

subplot(2,2,2); hold on; box on
plot(z*1e6,c2a(idx_t,:),'LineWidth',1.5)
xlabel('z (\mum)'); ylabel('[c2] (mM)')
xlim([-12 12]); ylim([0 1.1*max(max(c2a))])

subplot(2,2,3); hold on; box on
plot(z*1e6,c12a(idx_t,:),'LineWidth',1.5)
xlabel('z (\mum)'); ylabel('[c12] (mM)')
xlim([-12 12]); ylim([0 1.1*max(max(c12a))])

subplot(2,2,4); hold on; box on
plot(zcoords*1e6,c12conv_line*1e9,'LineWidth',1.5)
xlabel('z (\mum)'); ylabel('[c12]conv (nmol/m^2)')
xlim([-12 12]); ylim([0 1.1*max(max(c12conv_line*1e9))])

