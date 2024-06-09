% 1D DD simulation   
% By Deokhwa Seo, 2024/05/1
clc; clear;

% Initialize parameters
q = 1.602e-19; %[C] Elementary charge
kb = 1.380662e-23; %[J/K] Boltzmann constnat
eps0 = 8.854e-12; %[F/m] Vacuum permittivity
eps_s = 11.7 * eps0; %[F/m] Silicon dielectric constant 
T = 300; %[K] Temperature 
V_thermal = kb * T / q; %[V] Thermal voltage  
ni = 1e10 * 1e6; %[m^-3] Silicon intrinsic carrier concentration
Npoints = 501; % Number of points in space
Nsd_list = [1e19 2.5e19 5e19 7.5e19 1e20]; %[cm^-3] source/drain doping
Nch_list = [1e16 2.5e16 5e16 7.5e16 1e17 2.5e17 5e17 7.5e17 1e18 2.5e18 5e18 7.5e18]; %[cm^-3] channel doping
Lsd = 20; %[nm] source/drain length
Lch_list = 30;  %[nm] channel length (for single channel length)
% Lch_list = 20:10:100; %[nm] channel length (for multiple channel lengths)

% Generate parameter set
[Nsd_grid, Nch_grid, Lch_grid] = ndgrid(Nsd_list, Nch_list, Lch_list);
paramset_list = [Nsd_grid(:), Nch_grid(:), Lch_grid(:)];

% Loop through parameters 
Nparam = size(paramset_list, 1);
tic;
parfor param_idx = 1:Nparam

    % Open output file
    file_name = sprintf("Nsd_%.1e-Nch_%.1e-Lch_%d.dat", paramset_list(param_idx, :));
    file = fopen(file_name , "w");
    fprintf(file, "# Vd, Id, potential, charge\n");
  
    Nd1 = paramset_list(param_idx, 1) * 1e6; % S/D doping convert to m^-3
    Nd2 = paramset_list(param_idx, 2) * 1e6; % Channel doping convert to m^-3
    Lch = paramset_list(param_idx, 3);
    
    % Define x, dx, and interface boundaries 
    x = linspace(0, 2 * Lsd + Lch, Npoints);
    dx = (x(2)-x(1)) * 1e-9;
    N = length(x);
    idx_interface = round([Lsd, Lsd + Lch] / (x(2)-x(1)));
    idx_region1 = [1:idx_interface(1), idx_interface(2):N];
    idx_region2 = idx_interface(1)+1:idx_interface(2)-1;
    N_iter = 100;
    
    % Newton-Raphson to solve the poisson equation to get the initial
    % potential profile 
    J = zeros(N, N, N_iter);
    r = zeros(N, N_iter);
    phi = zeros(N, N_iter+1);
    dphi = zeros(N, N_iter);
    
    % Initial value
    phi(idx_region1,1) = V_thermal * log(Nd1 / ni);
    phi(idx_region2,1) = V_thermal * log(Nd2 / ni);
    idx_region1 = idx_region1(2:end-1); % Exclude endpoints
    
    for i = 1:N_iter
        r([1, N], i) = phi([1, N], i) - V_thermal * log(Nd1/ni);
        r(idx_region1, i) = phi(idx_region1+1,i) - 2 * phi(idx_region1, i) + phi(idx_region1-1, i) + q*dx^2/eps_s .* Nd1 - q*dx^2/eps_s * ni .* exp(phi(idx_region1,i)./V_thermal);
        r(idx_region2, i) = phi(idx_region2+1,i) - 2 * phi(idx_region2, i) + phi(idx_region2-1, i) + q*dx^2/eps_s .* Nd2 - q*dx^2/eps_s * ni .* exp(phi(idx_region2,i)./V_thermal);
        
        J_diag = zeros(N, 1);
        J_diag_u = [0; ones(N-2, 1)];
        J_diag_l = [ones(N-2, 1); 0];
        J_diag([1, N], 1) = 1;
        J_diag(idx_region1) = -2 - q*dx^2/eps_s / V_thermal * ni .* exp(phi(idx_region1, i)./V_thermal);
        J_diag(idx_region2) = -2 - q*dx^2/eps_s / V_thermal * ni .* exp(phi(idx_region2, i)./V_thermal);
        J(:,:,i) = diag(J_diag) + diag(J_diag_u, 1) + diag(J_diag_l, -1);
        
        dphi(:,i) = - J(:,:,i) \ r(:,i);
        phi(:, i+1) = phi(:, i) + dphi(:,i);
    
       
        dphi_rel = max(abs(dphi(:,i)./phi(:, i)));

        if max(abs(dphi_rel)) < 1e-9
            End_idx = i + 1;
            break;
        end
    end
    
    phi_poisson = phi(:,End_idx);
    
    
    % Use initial phi calcualted above and solve the continuity equation to
    % obtain the intial charge profile
    J = zeros(N, N, N_iter);
    r = zeros(N, N_iter);
    phi = phi_poisson .* ones(N, N_iter+1);
    n = zeros(N, N_iter+1);
    dn = zeros(N, N_iter);
    n(:,1) = ni .* exp(phi_poisson ./ V_thermal);
    
    for i = 1:N_iter
        r([1, N], i) = n([1, N], i) - Nd1;
        r(2:N-1, i) = (n((2:N-1)+1,i) + n((2:N-1),i)) .* (phi((2:N-1)+1, i) - phi((2:N-1),i)) - 2 * V_thermal .* (n((2:N-1)+1, i)-n((2:N-1),i)) ...
                           - (n((2:N-1),i) + n((2:N-1)-1,i)) .* (phi((2:N-1), i) - phi((2:N-1)-1,i)) + 2 * V_thermal .* (n((2:N-1), i)-n((2:N-1)-1,i));
        
        J_diag = zeros(N, 1);
        J_diag_u = zeros(N-1,1);
        J_diag_l = zeros(N-1,1);
        J_diag([1, N], 1) = 1;
        J_diag(2:N-1) = (phi((2:N-1)+1, i) - phi((2:N-1),i)) - (phi((2:N-1), i) - phi((2:N-1)-1,i)) + 4 * V_thermal;
        J_diag_u(2:end) =  (phi((2:N-1)+1, i) - phi((2:N-1),i)) - 2 * V_thermal;
        J_diag_l(1:end-1) =  - (phi((2:N-1), i) - phi((2:N-1)-1,i)) + 2 * V_thermal;
        J(:,:,i) = diag(J_diag) + diag(J_diag_u, 1) + diag(J_diag_l, -1);
        
        dn(:,i) = - J(:,:,i) \ r(:,i);
        n(:, i+1) = n(:, i) + dn(:,i);
    
        
        dn_rel = max(abs(dn(:,i)./n(:, i)));

        if dn_rel < 1e-9
            End_idx = i +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv1;
            break;
        end
    end
    
    n_continuity = n(:,End_idx);
    
    
    % Solve Poission and continuity equations simultaneously to using the initial 
    % potential and charge guesses calculated previously
    
    Vd = 0:0.05:1.2; % Simulating applied voltage range
    I = zeros(length(Vd),1);

    for Vd_idx = 1:length(Vd)
        J = zeros(2*N,  2*N, N_iter);
        r = zeros(2*N, N_iter);
        phi = zeros(N, N_iter+1);
        n = zeros(N, N_iter+1);
        phi(:,1) = phi_poisson;
        n(:,1) = n_continuity;
        idx_region1_poisson = 2*idx_region1 - 1;
        idx_region2_poisson = 2*idx_region2 - 1;
    
        for i = 1:N_iter
            % r for poisson
            r(2*[1, N]-1, i) = phi([1, N], i) - V_thermal * log(Nd1/ni) - [0; Vd(Vd_idx)];
            r(idx_region1_poisson, i) = phi(idx_region1+1, i) - 2 * phi(idx_region1, i) + phi(idx_region1-1, i) + q*dx^2/eps_s .* Nd1 - q*dx^2/eps_s * n(idx_region1, i);
            r(idx_region2_poisson, i) = phi(idx_region2+1, i) - 2 * phi(idx_region2, i) + phi(idx_region2-1, i) + q*dx^2/eps_s .* Nd2 - q*dx^2/eps_s * n(idx_region2, i);
    
            %r for continuity 
            idx_cont = 2:N-1;
            r(2 * [1, N], i) = n([1, N], i) - Nd1;
            r(2 * idx_cont , i) = n(idx_cont+1, i).* B((phi(idx_cont+1, i) - phi(idx_cont, i))/V_thermal) ...
                                  - n(idx_cont, i).* B((phi(idx_cont, i) - phi(idx_cont+1, i))/V_thermal) ...
                                  - n(idx_cont, i).* B((phi(idx_cont, i) - phi(idx_cont-1, i))/V_thermal) ...
                                  + n(idx_cont-1,i).* B((phi(idx_cont-1, i) - phi(idx_cont, i))/V_thermal);                        
    
            % J for poisson
            J(1, 1, i) = 1;
            J(2*N-1, 2*N-1, i) = 1;
            for j = 2:N-1
                J(2*j-1, (2*j-3):(2*j+1), i) = [1, 0, -2, - q*dx^2/eps_s, 1];
            end
    
            % J for continuity 
            J(2, 2, i) = 1;
            J(2*N, 2*N, i) = 1;
            for j = 2:N-1
                J(2*j, (2*j-3):(2*j+2), i) = [n(j, i) .* dB((phi(j,i) - phi(j-1, i))/V_thermal)/V_thermal + n(j-1, i) .* dB((phi(j-1,i) - phi(j, i))/V_thermal)/V_thermal; %phi(i-1)
                                          B((phi(j-1,i) - phi(j, i))/V_thermal); %n(i-1)
                                          - n(j+1, i) .* dB((phi(j+1,i) - phi(j, i))/V_thermal)/V_thermal - n(j, i) .* dB((phi(j,i) - phi(j+1, i))/V_thermal)/V_thermal...
                                          - n(j, i) .* dB((phi(j,i) - phi(j-1, i))/V_thermal)/V_thermal - n(j-1, i) .* dB((phi(j-1,i) - phi(j, i))/V_thermal)/V_thermal; %phi(i)
                                          - B((phi(j,i) - phi(j+1, i))/V_thermal) - B((phi(j,i) - phi(j-1, i))/V_thermal); %n(i)
                                          n(j+1, i) .* dB((phi(j+1,i) - phi(j, i))/V_thermal)/V_thermal + n(j, i) .* dB((phi(j,i) - phi(j+1, i))/V_thermal)/V_thermal; %phi(i+1)
                                          B((phi(j+1,i) - phi(j, i))/V_thermal)]'; %n(i+1)
            end 
            
            % Scale matrix 
            Cvector = zeros(2*N, 1);
            Cvector(1:2:2*N-1,1) = V_thermal;
            Cvector(2:2:2*N, 1) = Nd1;
            Cmatrix = spdiags(Cvector, 0, 2*N, 2*N);
            J_scaled = J(:,:,i) * Cmatrix;
            Rvector = 1./sum(abs(J_scaled),2);
            Rmatrix = spdiags(Rvector, 0 , 2*N, 2*N);
            J_scaled = Rmatrix * J_scaled;
            r_scaled = Rmatrix * r(:,i);
            update_scaled = - J_scaled \ r_scaled;
            update_vector = Cmatrix * update_scaled;
            
            dphi = update_vector(1:2:2*N-1);
            dn = update_vector(2:2:2*N);
            phi(:, i+1) = phi(:,i) + dphi;
            n(:, i+1) = n(:, i) + dn;
            
            % Find converged norms
            dn_rel = max(abs(dn./n(:, i)));
            dphi_rel = max(abs(dphi./phi(:, i))); 
            
            if dn_rel < 1e-9 && dphi_rel < 1e-9
                End_idx = i + 1;
                break;
            end
        end

        
        % Calculate current at the drain endpoint
        Id = q * 1417 * ((n(N,End_idx)+n(N-1,End_idx))/2*(phi(N,End_idx)-phi(N-1,End_idx))/dx - V_thermal*(n(N,End_idx)-n(N-1,End_idx))/dx)/1e8;
      
        % Save results
        output_line = [Vd(Vd_idx), Id, phi(:, i+1)', n(:, i+1)'];        
        fprintf(file, "%.6e ", output_line);
        fprintf(file, "\n");
    end
    
    fclose(file);
end
time_elapsed = toc;


% Bernoulli function       
function B = B(x)
    B = x./(exp(x)-1);
    B(find(abs(x)<1e-10)) = 1;
end

% Derivative of Bernoulli function
function dB = dB(x)
    dB = ((1-x) .* exp(x) - 1)./(exp(x)-1).^2;      
    dB(find(abs(x)<1e-6)) = -0.5;
end