function z_array=measurementsGenerationDistance(nodes_p,R,edge_matrix,rand_seed)

%Computes the measurements for each node

%Author: Angel F. Garcia Fernandez


dim_z=1;
N=size(nodes_p,2);
chol_R=chol(R)';

z_array=zeros(N,N,dim_z);

%Set lower triangular elements of edge_matrix to zero.
tri_edge_matrix = triu(edge_matrix);
[index_i,index_j]=find(tri_edge_matrix);


for i=1:length(index_i)
   
    pos_i=nodes_p(:,index_i(i));
    pos_j=nodes_p(:,index_j(i));
    
    %Measurement
    d_ij=sqrt(sum((pos_i-pos_j).^2));
    
    z=d_ij+chol_R*randn(rand_seed,1,1);
    z_array(index_i(i),index_j(i),:)=z;
    z_array(index_j(i),index_i(i),:)=z; %I repeat elements. This is not very efficient but convenient.
    
 
end





