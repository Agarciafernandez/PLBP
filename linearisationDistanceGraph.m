function lin_parameters=linearisationDistanceGraph(x_u,P_u,edge_matrix,weights)

%Computes the Statistical Linear Regression Parameters A,b,Omega at the initial iteration (which
%implies that the joint distribution of two nodes is independent as
%considered in the problem formulation)

%Author: Angel F. Garcia Fernandez

lin_parameters=struct;
%Set lower triangular elements of edge_matrix to zero.
tri_edge_matrix = triu(edge_matrix);
[index_i,index_j]=find(tri_edge_matrix);

for i=1:length(index_i)
    x_u_i=x_u(:,index_i(i));
    P_u_i=P_u(:,:,index_i(i));
    x_u_j=x_u(:,index_j(i));
    P_u_j=P_u(:,:,index_j(i));
    
    %We assume they are independent
    x_tot=[x_u_i;x_u_j];
    P_tot=[P_u_i,zeros(2);zeros(2),P_u_j];
    [A,b,Omega]=SLR_DistanceGraph(x_tot,P_tot,weights);
    
    lin_parameters(index_i(i),index_j(i)).A=A;
    lin_parameters(index_i(i),index_j(i)).b=b;
    lin_parameters(index_i(i),index_j(i)).Omega=Omega;
    
end

