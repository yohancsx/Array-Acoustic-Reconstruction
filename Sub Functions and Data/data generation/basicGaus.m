function gaussianOutput = basicGaus(x,y,x0,y0,amp,sigx,sigy,vo)
%basically computes a 3d gaussian, but with the special addition
%that values greater than 360 go to negative 360 and values less than 360 go
%to positive, us a sawtooth wave for this
xInput = pi*sawtooth((1/2)*(x - x0) + pi/2,1/2);
yInput = pi*sawtooth((1/2)*(y - y0) + pi/2,1/2);
gaussianOutput = abs(amp*exp( -1*( (((xInput).^2)./(2*sigx^2)) + (((yInput).^2)./(2*sigy^2)) ) ) + vo);
end

