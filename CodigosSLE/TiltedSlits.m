% Programa principal para plotear la curva
function [X, Y] = TiltedSlits(N, kappa)
  h = 0;
  dt = 1;
  dk = 0;
  for i = 1:N %Se plotean N puntos
    dk = delta_normal(dt, kappa); % Funcion directora por dist gaussiana
    %dk = delta_randomwalk(dt,kappa); % Funcion directora por random walk
    h = [0; Mapeo(h, dk, kappa, dt)];
  end

  X = real(h);
  Y = imag(h); % Pasamos al plano real en dos dimensiones para la figura

  figure(1),
  plot(X,Y);
  axis equal
  grid on
end


