function [A_sd, B_sd_n, B_sd_p] = matrices_solidphase(param)
%%%%%% ODE matrix formulation for Solid Phase of all cells in module %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% ODE: State-space formulation
A_mat = zeros(param.Nr-1);
B_mat = zeros(param.Nr-1,1);
for k = 1:(param.Nr-1)
    if k == 1
        A_mat(k,k) = -2;
        A_mat(k,k+1) = 2;
        B_mat(k,1) = 0;
    elseif k == (param.Nr-1)
        A_mat(k,k-1) = 2;
        A_mat(k,k) = -2;
        B_mat(k,1) = 2*(1+1/k);
    else
        A_mat(k,k-1) = 1*(1-1/k);
        A_mat(k,k) = -2;
        A_mat(k,k+1) = 1*(1+1/k);
        B_mat(k,1) = 0;
    end
end

% Creating a block-diagonal matrix with A matrices for all cells

%kron replaces each non-zero element of Nc x Nc identity matrix w/ (Nr - 1)x(Nr - 1) A_mat
A_sd = sparse(kron(eye(param.Nc),A_mat)); 

% Creating a column vector with B vector for all cells
%Repmat tiles Nc x 1 vector w/ (Nr-1)x 1 vector of B matrices for each cell
beta_n = 1./(param.F*param.A.*param.Ln.*param.a_sn.*param.delta_xn);
beta_p = 1./(param.F*param.A.*param.Lp.*param.a_sp.*param.delta_xp);

% Update for pack-level cell-to-cell variation:

B_sd_n = [];
B_sd_p = [];

for i = 1:param.Nc
    n_dummy = -beta_n(i)*B_mat;
    B_sd_n = [B_sd_n; n_dummy];
    
    p_dummy = beta_p(i)*B_mat;
    B_sd_p = [B_sd_p; p_dummy];
end

B_sd_n = sparse(B_sd_n);
B_sd_p = sparse(B_sd_p);

end
