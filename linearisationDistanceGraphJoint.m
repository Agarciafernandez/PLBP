function lin_parameters=linearisationDistanceGraphJoint(x_u_joint,P_u_joint,edge_matrix,pesos)

%Computes the Statistical Linear Regression Parameters A,b,Omega without
%assuming that the nodes are independent (which only happens at the first
%iteration)

%Author: Angel F. Garcia Fernandez


lin_parameters=struct;
%Set lower triangular elements of edge_matrix to zero.
tri_edge_matrix = triu(edge_matrix);
[index_i,index_j]=find(tri_edge_matrix);

for i=1:length(index_i)
    x_tot=x_u_joint(:,i);
    P_tot=P_u_joint(:,:,i);
    
    [A,b,Omega]=SLR_DistanceGraph(x_tot,P_tot,pesos);
    
    lin_parameters(index_i(i),index_j(i)).A=A;
    lin_parameters(index_i(i),index_j(i)).b=b;
    lin_parameters(index_i(i),index_j(i)).Omega=Omega;
    
end

