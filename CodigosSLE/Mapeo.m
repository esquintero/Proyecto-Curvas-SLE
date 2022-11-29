% Aqui se hace la funcion de mapeo para obtener los complejos
function hz = Mapeo(z, delta, kappa, dt)
  v = delta*delta/dt;
  alpha = 0;
  if delta > 0
    alpha = 0.5 - 0.5*sqrt(v / (16 + v));
  else
    alpha = 0.5 + 0.5*sqrt(v / (16 + v));
  end
  A = (z + 2 * sqrt(dt * (1 - alpha) / alpha)) .^ (1 - alpha);
  B = (z - 2 * sqrt(dt * alpha / (1 - alpha))) .^ alpha;
  hz = A .* B;
end

