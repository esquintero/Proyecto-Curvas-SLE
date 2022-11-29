% Obtencion de la funcion directora por distribucion gaussiana
function dk = delta_normal(dt,kappa)
  mu = 0; % Media de la distribucion
  sigma = sqrt(kappa*dt); % Desviacion estandar de la distribucion
  dk = sigma*randn(1)+mu;
