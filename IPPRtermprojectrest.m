%%%% Test ideas for IPPR term project
clear,clc;

% [X,Y] = meshgrid(-2:0.3:2,-2:0.3:2);
% ZX = cos(X) .* cos(Y);
% 
% a(1) = axes('position',[0.1 0.1 0.2 0.2]);
% a(2) = axes('position',[0.8 0.1 0.2 0.2]);
% a(3) = axes('position',[0.8 0.8 0.2 0.2]);
% a(4) = axes('position',[0.1 0.8 0.2 0.2]);
% a(5) = axes('position',[0.3 0.3 0.5 0.5]);
% 
% for ii = 1:5
%     axes(a(ii));
%     surf(ZX);
%     if ii == 1
%         view(37.5,30);
%     elseif ii == 2
%         view(-37.5,70);
%     elseif ii == 3
%         view(10, 30);
%     elseif ii == 4
%         view(0,20);
%     end
% end

% I = imread('IMG_2271.jpg');
% imshow(I,[]);
%caxis([20,60]);
% 
% figure(1);
% sphere;
% sp = findobj(1,'type','surface');
% cb1 = ['set(sp,''linestyle'',''none'')'];
% cb2 = ['set(sp,''linestyle'',''--'')'];
% cb3 = ['set(sp,''linestyle'','':'')'];
% cb4 = ['set(sp,''linestyle'',''-'')'];
% %%%% define context menu
% cmenu = uicontextmenu;
% set(sp,'uicontextmenu',cmenu);
% menp = uimenu(cmenu,'label','linetypes');
% item1 = uimenu(menp,'label','none','callback',cb1);
% item2 = uimenu(menp,'label','dashed','callback',cb2);
% item3 = uimenu(menp,'label','dotted','callback',cb3);
% item4 = uimenu(menp,'label','solid','callback',cb4);

% figure(1);
% sphere(25);
% x = [-2 -2 2 2];
% y = [-2 2 2 -2];
% z = [-2 -2 -2 -2];
% c = [-2 0 1 2];
% p1 = patch(x,y,z,c,'facelighting','none');
% xlabel('x');
% set(gca,'visible','off');

% [x,y,z] = sphere;
% ss = surface(x,y,z);
% view(3);
% set(ss,'facecolor','texturemap');


% x1 = [1 1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 -1 -1];
% y1 = [1 -1 -1 1 1 1 1 1 1 -1 -1 1 -1 -1 -1 -1];
% z1 = [1 1 1 1 1 -1 -1 1 -1 -1 -1 -1 -1 1 1 -1];
% clf;
% p = plot3(x1,y1,z1);
% set(p,'linewidth',3,'color','b');
% [XX,YY,ZZ] = sphere(15);
% hold on;
% h1 = mesh(XX,YY,ZZ);
% set(h1,'edgecolor','b','facecolor','c');
% h2 = mesh(2.*XX,2.*YY,2.*ZZ);
% set(h2,'edgecolor','r','facecolor','none');
% set(gca,'visible','off');
% axis square;
% view(3);

% figure;
% plot(x1,y1);view(3);

% figure(1); clf;
% subplot(2,1,1);
% surf(peaks);
% axis('off');
% title('unlit surface');
% 
% subplot(2,1,2);
% sp = surf(peaks);
% set(sp,'facecolor','red');
% set(sp,'edgecolor','none');
% set(sp,'facelighting','phong');
% light('position',[0 0 10],'style','local');
% %axis('off');
% title('lit surface');


% msgbox('hello','HELLO','warn','nonmodal');

%%%% create movie picture
% 
% x = -8:0.5:8;
% [XX,YY] = meshgrid(x);
% r = sqrt(XX.^2 + YY.^2) + eps;
% Z = sin(r) ./ r;
% surf(Z);
% theAxes = axis;
% %fmat = moviein(20);
% for ii = 1:20
%     surf(sin(2*pi*ii/20)*Z,Z);
%     axis(theAxes);
%     axis('off');
%     %fmat(:,ii) = getframe;
%     fmat(ii) = getframe;
% end
% movie(fmat,10);

x = linspace(0,10000,10000);
y = sin(x).*cos(x);
sound(y);



