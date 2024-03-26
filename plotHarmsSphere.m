function [H_Harms,V_Harms,handle_array]=plotHarmsSphere(N,radius,drawplot,formatstr)
%plot a Harms sphere's great circles, N-1 circles (so N segments) per dimension, return handles
H=linspace(0,pi*(N-1)/N,(N-1)); % don't need to do both 0 and pi since same circle
V=linspace(0,pi*(N-1)/N,(N-1));
handles = gobjects(N,2);
ii=1;
for azH=H 
  [~, handles(ii,1)] = plotCircle3D([0 0 0],[-sin(H) cos(H) 0],1,1,'k'); %ignore circlepts output
end
ii=1;
for elV=V % H and V arrays are same for Harms, but keeping this breakout for consistency with code for Fick
  [~, handles(ii,1)] = plotCircle3D([0 0 0],[-sin(V) 0 cos(V)],1,1,'b'); %ignore circlepts output
end
if nargout >= 1; H_Harms = H; end
if nargout >= 2; V_Harms = V; end
if nargout >= 3; handle_array = handles; end

end