function y = cos_ramp(x, fs, t_ramp)
% raise cosine ramp
% t_ramp: ramp time in s
up = (1+cos(pi+pi*linspace(0,1,round(t_ramp*fs))'))/2;
down = flipud(up); 
y = x;
y(1:round(fs*t_ramp)) = y(1:round(fs*t_ramp)).*up;
y(end-round(fs*t_ramp)+1:end) = y(end-round(fs*t_ramp)+1:end).*down; 
end