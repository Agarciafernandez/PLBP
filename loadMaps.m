function [pos_nodes,Nnodes,edge_matrix,x_p,P_p]=loadMaps(networkIndex,P_prior,P_prior_anchor,rand_seed)

%Author: Angel F. Garcia Fernandez
%Output
%pos_nodes: true positions of the nodes
%Nnodes: number of nodes
%Edge_matrix: represents the connectivity
%x_p: Prior means for each node
%P_p: Prior covariances for each node




filename1 = ['Maps/' num2str(networkIndex)];
load(filename1);
Nnodes=size(ZTOA,2);
pos_nodes=[Xagents;Xanchors]';
edge_matrix=double(ZTOA~=0);

%For drawing uncomment
% figure(1)
% plotNodes(nodes,edge_matrix)
% hold on
% plot(Xanchors(:,1),Xanchors(:,2),'xr','Linewidth',2,'MarkerSize',10)
% 
% hold off
% 
% axis([-10 110 0 100])
% axis equal
% grid on
% xlabel('x (m)')
% ylabel('y (m)')



x_p=zeros(2,Nnodes);
P_p=zeros(2,2,Nnodes);
for j=1:Nnodes
    
    if(j>size(Xagents,1))
        P_p(:,:,j)= P_prior_anchor;
        x_p(:,j)=pos_nodes(:,j)+chol(P_prior_anchor)'*randn(rand_seed,2,1);
        
    else
        P_p(:,:,j)= P_prior;
        x_p(:,j)=pos_nodes(:,j)+chol(P_prior)'*randn(rand_seed,2,1);
    end
end

