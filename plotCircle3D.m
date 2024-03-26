function [circlepts, handle] = plotCircle3D(center,normal,radius,drawplot,formatstr)
% plotCircle3D: Create a 100 point circle about centered on "center", about axis "normal", with radius radius
% center and normal can be either 3x1 or 1x3
% plots and returns handle if drawplot~=0, using formatstr like 'r:', otherwise just returns
% points

theta=(row(linspace(0,2*pi)))'; %force it to be a 100x1 vector (linspace defaults to 1x100)
v=null(row(normal)); %%force it to be a 1x3 vector; null() demands a row vector input, not col vector
center=row(center)'; %force it to be a 3x1 vector
circlepts=repmat(center,1,size(theta,1))+radius*(v(:,1)*cos(theta')+v(:,2)*sin(theta'));
if drawplot
  handle=plot3(circlepts(1,:),circlepts(2,:),circlepts(3,:),formatstr);
else
  handle=[];
end %if
end