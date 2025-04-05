% Jungho kim
% Spectral acceleration

function [PSA, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, gacc, dt)

gacc = gacc.*9.81.*1000;
% vel = cumtrapz(gacc)*dt;
% disp = cumtrapz(vel)*dt;
% Spectral solution
for i = 1:length(sPeriod)
    omegan = 2*pi/sPeriod(i);
    C = 2*xi*omegan;
    K = omegan^2;
    y(:,1) = [0;0];
    A = [0 1; -K -C];
    Ae = expm(A*dt);
    AeB = A\(Ae-eye(2))*[0;1];
    
    for k = 2:numel(gacc)
        y(:,k) = Ae*y(:,k-1) + AeB*gacc(k);
    end
    
    displ = (y(1,:))'./1000;                        
    veloc = (y(2,:))'./1000;                        
    foverm = omegan^2*displ;                        
    absacc = -2*xi*omegan*veloc - foverm;           
    
    % Extract peak values
    displ_max(i) = max(abs(displ));                 % Spectral relative displacement (m)
    veloc_max(i) = max(abs(veloc));                 % Spectral relative velocity (m/s)
    absacc_max(i) = max(abs(absacc));               % Spectral absolute acceleration (m/s2)
    
    foverm_max(i) = max(abs(foverm));               % Spectral value of lateral resisting force over mass (m/s2)
    pseudo_acc_max(i) = displ_max(i)*omegan^2;      % Pseudo spectral acceleration (m/s2)
    pseudo_veloc_max(i) = displ_max(i)*omegan;      % Pseudo spectral velocity (m/s)
    
    PSA(i) = pseudo_acc_max(i);                     % PSA (m/s2)
    PSV(i) = pseudo_veloc_max(i);                   % PSV (m/s)
    SA(i)  = absacc_max(i);                         % SA (m/s2)
    SV(i)  = veloc_max(i);                          % SV (m/s)
    SD(i)  = displ_max(i);                          % SD  (m)
  
    % Time series of acceleration, velocity and displacement response of
    % SDF oscillator
    OUT.acc(:,i) = absacc;
    OUT.vel(:,i) = veloc;
    OUT.disp(:,i) = displ;
end % function end


