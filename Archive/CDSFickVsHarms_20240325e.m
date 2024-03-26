% CDSFickVsHarms_YYYYMMDDx.m 
% Functions for plotting Fick and Harms spheres and wall projections, then
% computing errors that would occur if you use a Harms projection when
% calibrating a VOG system that you think is reported Fick H,V,T
%
% by Charley Della Santina 2024-03-xx
% R:\VNEL Common Files\Matlab VNEL Library\CDS Matlab routines\CDSFickVsHarms
% Still to do:
% 1. Make pretty, add labels, etc
% 2. Sensitivity analysis for offset of eye center relative to origin
% 3. Sensitivity analysis for offset of eye center relative to origin
% 4. Model "false torsion" effects

%{
References:
1) http://schorlab.berkeley.edu/passpro/oculomotor/html/chapter_3.html
2) Haslwanter T: Mathematics of three-dimensional eye rotations. Vision Res 35: 1727-1739, 1995
    https://jlc.jst.go.jp/DN/JALC/00106152399?type=list&lang=en&from=J-STAGE&dispptn=1
3) Takeshi Tsutsumi, Kei Tsukui, Ken Kitamura (2008) Translation of a coordinate system from Fick's frame of reference to rotational axes
    https://www.jstage.jst.go.jp/article/jser/67/6/67_6_522/_article
4) Moore et al (1994) Measurement of three dimensional eye position using image processing: a geometric approach
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=413351
5)
%}

NfacesFick=36;  % Make Fick sphere with 36x36 faces, so every 10 deg in azimuth and every 5 deg in Fick elevation
NHarms=18;      % Make Harms sphere with great circle every 10 deg
radius=1;   % 1 meter diameter of sphere. This isn't the eye's radius. It's the radius of a Fick/Harms sphere.
D=2;        % 2 meters from center of sphere to projection wall
makenewfig=1;   %tells plotFickSphereAndProjection() to make a new figure rather than add to existing

% compute and plot points on Fick sphere and Fick projection screen
% Note that for NfaceFick =36, H_Fick steps every 10 deg but V_Fick steps every 5 deg, not every 10 
[H_Fick,V_Fick,XsF,YsF,ZsF,XpF,YpF,ZpF,handles]=plotFickSphereAndProjection(NfacesFick,radius,D,makenewfig);

% compute and plot points on Harms sphere and Harms projection screen. 
% Note that for NHarms =18, H_Harms steps every 10 deg and V_Harms steps every 10 deg 
[H_Harms,V_Harms,XsH,YsH,ZsH,XpH,YpH,ZpH,handles]=plotHarmsSphere(NHarms,radius,D,makenewfig,'b:');

% Make a combined figure on same axes
plotFickSphereAndProjection(NfacesFick,radius,D,makenewfig);
plotHarmsSphere(NHarms,radius,D,~makenewfig,'b:');
title('Fick (black) and Harms (blue)');

% Make a second combined figure using the combined function made before I separated
% plotFickSphereAndProjection() and plotHarmsSphere()
overlayFickAndHarms();

% Compute, tabulate and plot the H_Fick and V_Fick you'd get if you used a
% Harms projection grid rather than a Fick projection grid for calibraton,
% mistakenly thinking you were using a Fick projection grid.
[H_Fick_ifHarmsgrid,V_Fick_ifHarmsgrid,~]=cart2sph(XpH,YpH,ZpH); %Matlab's standard function computes Fick azimuth and elevation (BUG: elevation SHOULD REALLY BE -V_Fick)

function [AZ,EL,Xs,Ys,Zs,Xp,Yp,Zp,handle]=plotFickSphereAndProjection(Nfaces,radius,D,makenewfig)
% [AZ,EL,Xs,Ys,Zs,Xp,Yp,Zp,handle]=plotFickSphereAndProjection(Nfaces,radius,D,makenewfig)
% 
% Plot an Nface Fick sphere and corresponding wall projection at x=D (typically D=2)
% AZ, EL are Fick angles H and V
% Xs,Y,Z are the intersection points in cartesian coordinates, which cam be
% converted to Fick using Matlab standard function cart2sph()
%======================================================= 
% Plot Fick sphere and rays to projection screen at x=D
% color red for az OR el >30deg
%======================================================= 
arguments
  Nfaces=36;  %number of faces on sphere; 36 is good for Fick, 18 is ok
  radius=1;   %radius of plotted sphere
  D=2;        %distance along x axis from origin to yz projection screen
  makenewfig=1; %1=create new figure, 0=plot into preexisting figure/axes
end  
  
Xslice=0.75; %will only show rays for part of sphere with Xs>Xslice
rad30deg=30*pi/180; % 30 degrees conerted to radians
rad20deg=20*pi/180; % 20 degrees conerted to radians
rad10deg=10*pi/180; % 10 degrees conerted to radians
maxangle=25; %degrees

%plot a Fick sphere of radius R and Nfaces faces
if makenewfig
    f1=figure;
end
[Xs, Ys, Zs] = sphere(Nfaces);
sphere(Nfaces);
handle=gca;
axis equal;
ggg=gca;
set(ggg,'Colormap',[1 1 1]); %sets colormap to 1
set(ggg.Children,'EdgeColor',[0.5 0.5 0.5]); 
hold on;
%for each vertex on the sphere with Xs > Xslice, plot a line in 3D from
%sphere origin to a plane at x=D
[AZ,EL,~]=cart2sph(Xs,Ys,Zs);
roundAZdeg=round(AZ*180/pi);
roundELdeg=round(EL*180/pi);
indexes = find(Xs>Xslice);
%indexes = find((abs(roundAZdeg)<maxangle) & (abs(roundELdeg)<maxangle));

for ii=1:length(indexes)
  jj=indexes(ii);
  [az,el,r]=cart2sph(Xs(jj), Ys(jj), Zs(jj)); %az = +H_Fick, el= -V_Fick
  Xp=D; %p for projection screen
  Yp=D*tan(az);
  Zp=sqrt(Xp*Xp+Yp*Yp)*tan(el);

  color=[0 0 0]; %just make all the Fick black for now
  %change color to signify angular distance from primary position
%{
  eccentricity=acosd(Xs(jj)); %how far from primary position, in deg
  if eccentricity >= 40
    color=[1 0 0];
  elseif eccentricity <= 20
    color=[0 0.75 0];
  else
    color=[0 0 1];
  end
%}
  plot3([Xs(jj);Xp],[Ys(jj);Yp],[Zs(jj);Zp],'Color',color);
  plot3(Xp,Yp,Zp,'Marker','o','MarkerSize',3,'MarkerFaceColor',color,'Color',color);
end

%Make plot pretty
htitle=title('Fick');
hxlabel=xlabel('x');
hylabel=ylabel('y');
hzlabel=zlabel('z');
h1=plot3([1 2],[0 0],[0 0],'g*-','LineWidth',3);
h2=plot3([0 0],[1 2],[0 0],'b*-','LineWidth',3);
h3=plot3([0 0],[0 0],[1 2],'r*-','LineWidth',3);
hx=text(2.2,0,0,'x','FontSize',16);
hy=text(0,2.2,0,'y','FontSize',16);
hz=text(0,0,2.2,'z','FontSize',16);
haz20 =text(1.05*cos(20*pi/180),1.05*sin( 20*pi/180),0,'+20','FontSize',8);
haz20n=text(1.05*cos(20*pi/180),1.05*sin(-20*pi/180),0,'-20','FontSize',8);
hel20 =text(1.05*cos(20*pi/180),0,1.05*sin( 20*pi/180),'+20','FontSize',8);
hel20n=text(1.05*cos(20*pi/180),0,1.05*sin(-20*pi/180),'-20','FontSize',8);
haz40 =text(1.05*cos(40*pi/180),1.05*sin( 40*pi/180),0,'+40','FontSize',8);
haz40n=text(1.05*cos(40*pi/180),1.05*sin(-40*pi/180),0,'-40','FontSize',8);
hel40 =text(1.05*cos(40*pi/180),0,1.05*sin( 40*pi/180),'+40','FontSize',8);
hel40n=text(1.05*cos(40*pi/180),0,1.05*sin(-40*pi/180),'-40','FontSize',8);
axis equal;
axis vis3d;

end

%=======================================================================
function [H_Harms,V_Harms,XsH,YsH,ZsH,XpH,YpH,ZpH,handle_array]=plotHarmsSphere(N,radius,D,makenewfig,formatstr)
% [H_Harms,V_Harms,X,Y,Z,handle_array]=plotHarmsSphere(N,radius,drawplot,formatstr)
% 
% Plot a Harms sphere's great circles, N-1 circles (so N segments) per dimension, return handles
% H_Harms,V_Harms are the Harms_azimuth and Harms_elevation 
%   ***BUG: Positive V_Harms is in the top half of the sphere here; should probably invert it before computing rotation vectors  
% X,Y,Z are Harms great circle intersection points in cartesian coordinates, which can be converted to Fick using Matlab standard function cart2sph()
% 
arguments
  N=18;       % Draw N Harms great circles, so spacing is every 10 deg
  radius=1;   % radius of plotted sphere
  D=2;        % distance along x axis from origin to yz projection screen       
  makenewfig=1;     %plots the results it drawplot ~=0
  formatstr='b:'  %defaults appearance for Harms grid on sphere
end

if makenewfig
    f2=figure;
end

Xslice=0.75; %will only show rays for part of sphere with Xs>Xslice
rad30deg=30*pi/180; % 30 degrees conerted to radians
rad20deg=20*pi/180; % 20 degrees conerted to radians
rad10deg=10*pi/180; % 10 degrees conerted to radians
maxangle=25; %degrees
drawplot=1;

H=linspace(0,pi*(N-1)/N,(N-1)); % don't need to do both 0 and pi since same circle
V=linspace(0,pi*(N-1)/N,(N-1));
handles = gobjects(N,2);
ii=1;
step=pi/N;
last=(N-1)*pi/N;
for hh=0:step:last
  [~, handles(ii,1)] = plotCircle3D([0 0 0],[-sin(hh) cos(hh) 0],radius,drawplot,formatstr); %ignore circlepts output
  ii=ii+1;
  if drawplot
    hold on
  end
end
ii=1;
for vv=0:step:last % H and V arrays are same for Harms, but keeping this breakout for consistency with code for Fick
  [~, handles(ii,2)] = plotCircle3D([0 0 0],[-sin(vv) 0 cos(vv)],radius,drawplot,formatstr); %ignore circlepts output
  ii=ii+1;
     if drawplot
    hold on
  end
end

%plot a sphere inside the Harms circles to obscure the back side
if drawplot
  [x,y,z]=sphere;
  hSurface=surf(x*0.99,y*0.99,z*0.99);
  set(hSurface,'FaceColor',[1 1 1],'EdgeColor','none')
  % make plot axes pretty
  axis equal
  axis vis3d
  htitle=title('Harms');
  hxlabel=xlabel('x');
  hylabel=ylabel('y');
  hzlabel=zlabel('z');
  grid on;
end

%compute the intersections and output them analogous to Fick sphere outputs
%given by standard Matlab function [x,y,z]=sphere(N)
N_HarmsFaces=N*N*2;
N_HarmsVerticesWithDuplicates=N*N*2;
N_HarmsVerticesNoDuplicates=(N-1)*(N-1)*2+4;
%dat=zeros(N,N,2,3); %dimensions: NH, NV, 2 intersections per circle pair, 3 dimensions (x,y,z) per intersection
dat=zeros(N,N,3); %dimensions: NH, NV, 3 dimensions (x,y,z) per intersection
ii=1;
for hh=0:step:last
  jj=1;
  for vv=0:step:last
    %intersection=null([sin(hh) -cos(hh) 0; sin(vv) 0 -cos(vv)]);
    u=cross([sin(hh) -cos(hh) 0],[sin(vv) 0 -cos(vv)]);
    intersection=radius*u/norm(u);
    dat(ii,jj,:)=intersection;
    jj=jj+1;
  end
  ii=ii+1;
end

%figure,
Xh=round([dat(:,:,1) -dat(:,:,1)],5);
Yh=round([dat(:,:,2) -dat(:,:,2)],5);
Zh=round([dat(:,:,3) -dat(:,:,3)],5);

plot3(Xh,Yh,Zh,'ok');
hold on;
[x,y,z]=sphere;
  hSurface=surf(x*0.99,y*0.99,z*0.99);
  set(hSurface,'FaceColor',[1 1 1],'EdgeColor','none')
axis equal
axis vis3d

%now plot rays to a Harms projection screen
%for each vertex on the sphere with Xs > Xslice, plot a line in 3D from
%sphere origin to a plane at x=D
Xslice=0.75;
indexes = find(Xh>Xslice);
for ii=1:length(indexes)
  jj=indexes(ii);
  [az,el,r]=cart2sph(Xh(jj), Yh(jj), Zh(jj)); %az = +H_Fick, el= -V_Fick
  Xp=D; %p for projection screen
  Yp=D*tan(az);
  Zp=sqrt(Xp*Xp+Yp*Yp)*tan(el);
  
  color=[0 0 1]; %just make all the Harms blue for now
  %change color to signify angular distance from primary position
%{
  eccentricity=acosd(Xs(jj)); %how far from primary position, in deg
  if eccentricity >= 40
    color=[1 0 0];
  elseif eccentricity <= 20
    color=[0 0.75 0];
  else
    color=[0 0 1];
  end
%}
  plot3([Xh(jj);Xp],[Yh(jj);Yp],[Zh(jj);Zp],'Color',color);
  plot3(Xp,Yp,Zp,'Marker','o','MarkerSize',4,'MarkerFaceColor',color,'Color',color);
end

%Make plot pretty
htitle=title('Harms');
hxlabel=xlabel('x');
hylabel=ylabel('y');
hzlabel=zlabel('z');
h1=plot3([1 2],[0 0],[0 0],'g*-','LineWidth',3);
h2=plot3([0 0],[1 2],[0 0],'b*-','LineWidth',3);
h3=plot3([0 0],[0 0],[1 2],'r*-','LineWidth',3);
hx=text(2.2,0,0,'x','FontSize',16);
hy=text(0,2.2,0,'y','FontSize',16);
hz=text(0,0,2.2,'z','FontSize',16);
haz20 =text(1.05*cos(20*pi/180),1.05*sin( 20*pi/180),0,'+20','FontSize',8);
haz20n=text(1.05*cos(20*pi/180),1.05*sin(-20*pi/180),0,'-20','FontSize',8);
hel20 =text(1.05*cos(20*pi/180),0,1.05*sin( 20*pi/180),'+20','FontSize',8);
hel20n=text(1.05*cos(20*pi/180),0,1.05*sin(-20*pi/180),'-20','FontSize',8);
haz40 =text(1.05*cos(40*pi/180),1.05*sin( 40*pi/180),0,'+40','FontSize',8);
haz40n=text(1.05*cos(40*pi/180),1.05*sin(-40*pi/180),0,'-40','FontSize',8);
hel40 =text(1.05*cos(40*pi/180),0,1.05*sin( 40*pi/180),'+40','FontSize',8);
hel40n=text(1.05*cos(40*pi/180),0,1.05*sin(-40*pi/180),'-40','FontSize',8);
axis equal;
axis vis3d;

%set output variables
if nargout >= 1; H_Harms = H; end
if nargout >= 2; V_Harms = V; end
if nargout >= 3; XsH = Xh; end
if nargout >= 4; YsH = Yh; end
if nargout >= 5; ZsH = Zh; end
if nargout >= 6; XpH = Xp; end
if nargout >= 7; YpH = Yp; end
if nargout >= 8; ZpH = Zp; end
if nargout >= 9; handle_array = handles; end

end

function overlayFickAndHarms()
%======================================================= 
% Plot Fick sphere and rays to projection screen at x=D
% then overlay Harms sphere and rays.
%======================================================= 
Nfaces=36;  %number of faces on Fick sphere
N=18;       %number of great circles per dimension on Harms sphere
radius=1;   %radius of plotted sphere
D=2;        %distance along x axis from origin to yz projection screen
Xslice=0.75; %will only show rays for part of sphere with Xs>Xslice
drawplot=1; %draw the Harms plot
formatstr='b:'; %format string for harms plot

rad30deg=30*pi/180; % 30 degrees conerted to radians
rad20deg=20*pi/180; % 20 degrees conerted to radians
rad10deg=10*pi/180; % 10 degrees conerted to radians
maxangle=31; %degrees

%plot a Fick sphere of radius R and Nfaces faces
f1=figure;
[Xs, Ys, Zs] = sphere(Nfaces);
sphere(Nfaces);
axis equal;
set(gca,'Colormap',1+0.*(get(gca,'Colormap'))); %sets colormap to 1
hold on;
%for each vertex on the sphere with Xs > Xslice, plot a line in 3D from
%sphere origin to a plane at x=D
[AZ,EL,~]=cart2sph(Xs,Ys,Zs);
roundAZdeg=round(AZ*180/pi);
roundELdeg=round(EL*180/pi);
indexes = find(Xs>Xslice);
%indexes = find((abs(roundAZdeg)<maxangle) & (abs(roundELdeg)<maxangle));

for ii=1:length(indexes)
  jj=indexes(ii);
  [az,el,r]=cart2sph(Xs(jj), Ys(jj), Zs(jj)); %az = +H_Fick, el= -V_Fick
  Xp=D; %p for projection screen
  Yp=D*tan(az);
  Zp=sqrt(Xp*Xp+Yp*Yp)*tan(el);
  color=[0 0 0]+0.3;
%   if (abs(az) >= rad30deg) || (abs(el) >= rad30deg)
%     color=[1 0 0];
% %   elseif (abs(az) <= rad10deg) & (abs(el) <= rad10deg)
% %     color=[0 0.5 0];
% %   else
% %     color=[0 0 1];

  plot3([Xs(jj);Xp],[Ys(jj);Yp],[Zs(jj);Zp],'Color',color);
  plot3(Xp,Yp,Zp,'Marker','o','MarkerSize',3,'MarkerFaceColor',color,'Color',color);
end
axis vis3d;
htitle=title('Fick and Harms');
hxlabel=xlabel('x');
hylabel=ylabel('y');
hzlabel=zlabel('z');

H=linspace(0,pi*(N-1)/N,(N-1)); % don't need to do both 0 and pi since same circle
V=linspace(0,pi*(N-1)/N,(N-1));
handles = gobjects(N,2);
ii=1;
step=pi/N;
last=(N-1)*pi/N;
for hh=0:step:last
  [~, handles(ii,1)] = plotCircle3D([0 0 0],[-sin(hh) cos(hh) 0],radius,drawplot,formatstr); %ignore circlepts output
  ii=ii+1;
  if drawplot
    hold on
  end
end
ii=1;
for vv=0:step:last % H and V arrays are same for Harms, but keeping this breakout for consistency with code for Fick
  [~, handles(ii,2)] = plotCircle3D([0 0 0],[-sin(vv) 0 cos(vv)],radius,drawplot,formatstr); %ignore circlepts output
  ii=ii+1;
     if drawplot
    hold on
  end
end

%plot a sphere inside the Harms circles to obscure the back side
if drawplot
  [x,y,z]=sphere;
  hSurface=surf(x*0.99,y*0.99,z*0.99);
  set(hSurface,'FaceColor',[1 1 1],'EdgeColor','none')
  % make plot axes pretty
  axis equal
  axis vis3d
  htitle=title('Harms');
  hxlabel=xlabel('x');
  hylabel=ylabel('y');
  hzlabel=zlabel('z');
  grid on;
end

%compute the intersections and output them analogous to Fick sphere outputs
%given by standard Matlab function [x,y,z]=sphere(N)
N_HarmsFaces=N*N*2;
N_HarmsVerticesWithDuplicates=N*N*2;
N_HarmsVerticesNoDuplicates=(N-1)*(N-1)*2+4;
%dat=zeros(N,N,2,3); %dimensions: NH, NV, 2 intersections per circle pair, 3 dimensions (x,y,z) per intersection
dat=zeros(N,N,3); %dimensions: NH, NV, 3 dimensions (x,y,z) per intersection
ii=1;
for hh=0:step:last
  jj=1;
  for vv=0:step:last
    %intersection=null([sin(hh) -cos(hh) 0; sin(vv) 0 -cos(vv)]);
    u=cross([sin(hh) -cos(hh) 0],[sin(vv) 0 -cos(vv)]);
    intersection=radius*u/norm(u);
    dat(ii,jj,:)=intersection;
    jj=jj+1;
  end
  ii=ii+1;
end

%figure,
Xh=round([dat(:,:,1) -dat(:,:,1)],5);
Yh=round([dat(:,:,2) -dat(:,:,2)],5);
Zh=round([dat(:,:,3) -dat(:,:,3)],5);

plot3(Xh,Yh,Zh,'ok');
hold on;
[x,y,z]=sphere;
  hSurface=surf(x*0.99,y*0.99,z*0.99);
  set(hSurface,'FaceColor',[1 1 1],'EdgeColor','none')
axis equal
axis vis3d

%now plot rays to a Harms projection screen
%for each vertex on the sphere with Xs > Xslice, plot a line in 3D from
%sphere origin to a plane at x=D
Xslice=0.75;
indexes = find(Xh>Xslice);
for ii=1:length(indexes)
  jj=indexes(ii);
  [az,el,r]=cart2sph(Xh(jj), Yh(jj), Zh(jj)); %az = +H_Fick, el= -V_Fick
  Xp=D; %p for projection screen
  Yp=D*tan(az);
  Zp=sqrt(Xp*Xp+Yp*Yp)*tan(el);
  color=[0 0 1]; %just make all the Harms blue for now
  %change color to signify angular distance from primary position
%{
  eccentricity=acosd(Xh(jj)); %how far from primary position, in deg
  if eccentricity >= 40
    color=[1 0 0];
  elseif eccentricity <= 20
    color=[0 0.75 0];
  else
    color=[0 0 1];
  end
%}
  plot3([Xh(jj);Xp],[Yh(jj);Yp],[Zh(jj);Zp],'Color',color);
  plot3(Xp,Yp,Zp,'Marker','o','MarkerSize',4,'MarkerFaceColor',color,'Color',color);
end

%Make plot pretty
htitle=title('Fick (black) and Harms (blue)');
hxlabel=xlabel('x');
hylabel=ylabel('y');
hzlabel=zlabel('z');
h1=plot3([1 2],[0 0],[0 0],'g*-','LineWidth',3);
h2=plot3([0 0],[1 2],[0 0],'b*-','LineWidth',3);
h3=plot3([0 0],[0 0],[1 2],'r*-','LineWidth',3);
hx=text(2.2,0,0,'x','FontSize',16);
hy=text(0,2.2,0,'y','FontSize',16);
hz=text(0,0,2.2,'z','FontSize',16);
haz20 =text(1.05*cos(20*pi/180),1.05*sin( 20*pi/180),0,'+20','FontSize',8);
haz20n=text(1.05*cos(20*pi/180),1.05*sin(-20*pi/180),0,'-20','FontSize',8);
hel20 =text(1.05*cos(20*pi/180),0,1.05*sin( 20*pi/180),'+20','FontSize',8);
hel20n=text(1.05*cos(20*pi/180),0,1.05*sin(-20*pi/180),'-20','FontSize',8);
haz40 =text(1.05*cos(40*pi/180),1.05*sin( 40*pi/180),0,'+40','FontSize',8);
haz40n=text(1.05*cos(40*pi/180),1.05*sin(-40*pi/180),0,'-40','FontSize',8);
hel40 =text(1.05*cos(40*pi/180),0,1.05*sin( 40*pi/180),'+40','FontSize',8);
hel40n=text(1.05*cos(40*pi/180),0,1.05*sin(-40*pi/180),'-40','FontSize',8);
axis equal;
axis vis3d;

end

%=======================================================================
function v = col(v)
%Forces vector to be a column; complement to row(v), which forces a vector to be a row.
    v = row(v)';
end



%=======================================================================
function [circlepts, handle] = plotCircle3D(center,normal,radius,drawplot,formatstr)
% plotCircle3D: Create a circle about centered on "center", about axis "normal", with radius radius
% center and normal can be either 3x1 or 1x3
% plots and returns handle if drawplot~=0, using formatstr like 'b:', otherwise just returns points

theta=(row(linspace(0,2*pi)))'; %force it to be a 3x1 vector
v=null(row(normal)); %%force it to be a 1x3 vector; null() demands a row vector input, not col vector
center=row(center)'; %force it to be a 3x1 vector
circlepts=repmat(center,1,size(theta,1))+radius*(v(:,1)*cos(theta')+v(:,2)*sin(theta'));
if drawplot
  handle=plot3(circlepts(1,:),circlepts(2,:),circlepts(3,:),formatstr);
else
  handle=[];
end %if
end


%=======================================================================
function gridpts = make3DgridYZ(center,ystep,zstep,N)
% Make an N-point square grid array of points in 3D about given center in
% the ax1,ax2 plane with (oddrootN)^2 points spaced by step
 
% UNDER CONSTRUCTION
end

%=======================================================================
function [gridpts, vectorpts] = make3DgridARB(center,ax1,ax2,step,N)
% Make an N-point square grid array of points in 3D about given center in
% the ARBitrary ax1,ax2 plane with (sqrt(N))^2 points spaced by step,
% with ax1 and ax2 unitized vectors and step = number of unit steps along
% each per point spacing. ax1 and ax2 must be 3x1 vectors.

ax1=ax1/norm(ax1); %normalize to a unit vector
ax2=ax2/norm(ax2);

rootN=round(sqrt(N));
M=(rootN-1)/2;
%gridpts = mesh((-M:1:M),(M:1:-M)); 
pts1 = zeros(3,N); % all N points in a single 3xN array
pts2 = zeros(3,rootN,rootN); % the same points in a sqrt(N) x sqrt(N) array of 3x1 vectors
count=1
offset=-M*step*(ax1+ax2);
for row = 1:1:rootN
  for col = rootN:-1:1
    pts1(:,count)= ((row-1)*ax1 + (col-1)*ax2)*step + center + offset;
    count=count+1;
    pts2(:,row,col)= ((row-1)*ax1 + (col-1)*ax2)*step + center + offset;
  end
end
vectorpnts = [pts1(1,:)' pts1(2,:)' pts1(3,:)'];
gridpnts = pts2;
figure,plot3(vectorpnts(:,1),vectorpnts(:,2),vectorpnts(:,3),'r*');
end

%=======================================================================
function circlepts = make3Dcircle(center,ax,radius,N)
% Make an N-point circle in 3D about given center & axis, with radius

% First make an N-point unit circle in the xy plane about origin
t = linspace(0,2*pi);
x = cos(t);
y = sin(t);
z = 0*t;
pnts = [x;y;z];

% unit normal for original plane
n0 = [0;0;1]; 
n0 = n0/norm(n0);

% unit normal for plane to rotate into 
% plane is orthogonal to n1... given by equation
% n1(1)*x + n1(2)*y + n1(3)*z = 0
n1 = ax/norm(ax); 

% theta is the angle between normals
c = dot(n0,n1) / ( norm(n0)*norm(n1) ); % cos(theta)
s = sqrt(1-c*c);                        % sin(theta)
u = cross(n0,n1) / ( norm(n0)*norm(n1) ); % rotation axis...
u = u/norm(u); % ... as unit vector
C = 1-c;

% the rotation matrix
R = [u(1)^2*C+c,          u(1)*u(2)*C-u(3)*s,   u(1)*u(3)*C+u(2)*s
     u(2)*u(1)*C+u(3)*s,  u(2)^2*C+c,           u(2)*u(3)*C-u(1)*s
     u(3)*u(1)*C-u(2)*s,  u(3)*u(2)*C+u(1)*s,   u(3)^2*C+c        ];

% Rotated points
circlepts = R*pnts;
end

%{
================
== SCRAP TEXT ==
================

origin=[0;0;0]; %center of eye
R=1;            %m (meters), radius of eye sphere on which we will illustrate points
D=2;            %m, distance along +x axis from origin to projection screen perpendicular to +x axis
Reye=0.013      %m, eye radius from center of rotation to anterior surface is about 13mm
gridspacingJHOC6030=11*0.0254 %m, spacing between dots on grid is 11 inches
gridspacing = gridspacingJHOC6030;
laserhalfangledeg= 3.2 %deg, angle from Neurolign laser center dot to top, bottom or side of laser calibration cross
 
d_az = 10*pi/180;   %steps in azimuth, aka longitude, in radians
d_el = 10*pi/180;   %steps in azimuth, aka longitude, in radians


Built-in Matlab functions:

[x,y,z] = sph2cart(azimuth,elevation,r) % https://www.mathworks.com/help/matlab/ref/sph2cart.html
  x = r .* cos(elevation) .* cos(azimuth)
  y = r .* cos(elevation) .* sin(azimuth)
  z = r .* sin(elevation)

[azimuth,elevation,r] = cart2sph(x,y,z) % https://www.mathworks.com/help/matlab/ref/cart2sph.html
  azimuth = atan2(y,x)
  elevation = atan2(z,sqrt(x.^2 + y.^2))
  r = sqrt(x.^2 + y.^2 + z.^2)




%}