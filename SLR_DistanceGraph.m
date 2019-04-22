function   [A,b,Omega]=SLR_DistanceGraph(meank,Pk,weights)

%Computes SLR parameters A,b, Omega for a given edge which has the joint
%distribution over the nodes indicated by meank and Pk
%Author: Angel F. Garcia Fernandez

Nx=length(meank);
W0=weights(1);
Nz=1;

chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
sigma_points=repmat(meank,1,length(weights))+sigma_points;


sigma_points_z=sqrt((sigma_points(1,:)-sigma_points(3,:)).^2+(sigma_points(2,:)-sigma_points(4,:)).^2);

z_pred_ukf=sigma_points_z*weights';

var_pred_ukf=zeros(Nz);
var_xz_ukf=zeros(Nx,Nz);

for j=1:length(weights)
    sub_z_j=sigma_points_z(:,j)-z_pred_ukf;
    var_pred_ukf=var_pred_ukf+weights(j)*(sub_z_j*sub_z_j');
    var_xz_ukf=var_xz_ukf+weights(j)*(sigma_points(:,j)-meank)*sub_z_j';
end

%Statistical linearisaltion

A=var_xz_ukf'/Pk;
b=z_pred_ukf-A*meank;
Omega=var_pred_ukf-A*Pk*A';