%This code runs the  posterior linearisation belief propagation (PLBP) algorithm, based on sigma-points, proposed in

% A. F. García-Fernández, L. Svensson, S. Särkkä, "Cooperative localization using posterior linearization belief propagation",
% IEEE Transactions on Vehicular Technology, vol. 67, no. 1, pp. 832-836, Jan. 2018.


%Author: Angel F. Garcia Fernandez

clear
rand_s1 = RandStream.create('mrg32k3a','NumStreams',1,'seed',1);

%Prior covariance matrix (non-anchor nodes)
P_prior=diag([100,100]);
%Prior covariance matrix (anchor nodes)
P_prior_anchor=diag([0.01,0.01]);
%Measurement noise variance
R=1;
%Load scenario (obtained from Prof. Henk Wymeersch web page)
[pos_nodes,N,edge_matrix,x_p,P_p]=loadMaps(7,P_prior,P_prior_anchor,rand_s1);

%Nsteps: Number of posterior linearisation steps
Nsteps=5;
%Nsteps_loopy: Number of (loopy) belief propagation steps for each poterior
%linearisation iteration
Nsteps_loopy=5;

%Nmc: Number of Monte Carlo runs
Nmc=50;
%square_error_mc: Square error for each Monte Carlo run
square_error_mc=zeros(Nmc,Nsteps+1);


%Unscented Transform parameters
W0=1/3; %Weight of sigma-point at the mean
Nx=2*size(pos_nodes,1); %Dimension
Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)];

for i=1:Nmc
    
    tic
    rand_s1 = RandStream.create('mrg32k3a','NumStreams',1,'seed',i);
    
    
    %Generation of measurements
    z_array=measurementsGenerationDistance(pos_nodes,R,edge_matrix,rand_s1);
    x_u=x_p;
    P_u=P_p;
    
    square_error_mc(i,1)=sum(sum((pos_nodes-x_u).^2));
    
   
    
    for j=1:Nsteps
               
        %Linearisation        
        if(j==1)
            lin_parameters=linearisationDistanceGraph(x_u,P_u,edge_matrix,weights);
        else
            lin_parameters=linearisationDistanceGraphJoint(x_u_joint,P_u_joint,edge_matrix,weights);
        end
        %BP
        [x_u2,P_u2,x_u_joint,P_u_joint]=Loopy_BP_linearised_model_joint(lin_parameters,edge_matrix,x_p,P_p,z_array,R,Nsteps_loopy);       
        x_u=x_u2;
        P_u=P_u2;
        
        
        %     rms_nodes_it(j+1)=sqrt(sum(sum((pos_nodes-x_u).^2))/size(pos_nodes,2));
        square_error_mc(i,j+1)=sum(sum((pos_nodes-x_u).^2));
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' seg'])
    
end


rms_error_it=sqrt(sum(square_error_mc,1)/(size(pos_nodes,2)*Nmc));



figure(2)
plot(rms_error_it)
xlabel('Iterations')
ylabel('RMS (m)')
grid on
title('RMS posterior linearisation')


