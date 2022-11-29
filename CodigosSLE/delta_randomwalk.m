% Obtencion de la funcion directora por random walk
function dk = delta_randomwalk(dt, kappa)
  sign = rand(1);
  if sign > 0.5
    dk = sqrt(dt*kappa);
  else
    dk = -sqrt(dt*kappa);
  end
end
