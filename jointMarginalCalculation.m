function [x_joint,P_joint,x_u,P_u]=jointMarginalCalculation(messages,lin_parameters,z_array,edge_matrix,x_p,P_p,R)

%This function computes the joint distributions of the nodes that are
%connected as well as the marginal distribution of each node

%Author: Angel F. Garcia Fernandez


%Set lower triangular elements of edge_matrix to zero.
tri_edge_matrix = triu(edge_matrix);
[index_i,index_j]=find(tri_edge_matrix);

N=length(index_i);

x_joint=zeros(4,N);
P_joint=zeros(4,4,N);


x_u=x_p;
P_u=P_p;

marginal_calc=zeros(1,size(x_p,2));

for i=1:N
    
    node_i=index_i(i);
    node_j=index_j(i);
    
    %First we do node i (with its neighbours except node_j
    nodes_input=edge_matrix(:,node_i);
    nodes_input(node_j)=0; %We remove node_j
    list_nodes_input=find(nodes_input==1);
    %We compute the KF updates
    x_p_i=x_p(:,node_i);
    P_p_i=P_p(:,:,node_i);
    
    x_u_i=x_p_i;
    P_u_i=P_p_i;
    
    for p=1:length(list_nodes_input)
        node_p=list_nodes_input(p); %node_p: neighbour of transmitting node
        alfa_p=messages(node_p,node_i).alfa;
        Hmess_p=messages(node_p,node_i).Hmess;
        gamma_p=messages(node_p,node_i).gamma;
        [x_u_i,P_u_i]=linear_kf_update(x_u_i,P_u_i,Hmess_p,0,0,gamma_p,alfa_p);       
    end
    
    %We do the same with node j
    nodes_input=edge_matrix(:,node_j);
    nodes_input(node_i)=0; %We remove node_i
    list_nodes_input=find(nodes_input==1);
    %We compute the KF updates
    x_p_j=x_p(:,node_j);
    P_p_j=P_p(:,:,node_j);
    
    x_u_j=x_p_j;
    P_u_j=P_p_j;
    
    for p=1:length(list_nodes_input)
        node_p=list_nodes_input(p); %node_p: neighbour of transmitting node
        alfa_p=messages(node_p,node_j).alfa;
        Hmess_p=messages(node_p,node_j).Hmess;
        gamma_p=messages(node_p,node_j).gamma;
        [x_u_j,P_u_j]=linear_kf_update(x_u_j,P_u_j,Hmess_p,0,0,gamma_p,alfa_p);       
    end
    
    %Final update is (we need the linearisation)
    A=lin_parameters(node_i,node_j).A;
    b=lin_parameters(node_i,node_j).b;
    Omega=lin_parameters(node_i,node_j).Omega;
    z=squeeze(z_array(node_i,node_j,:));
    
    joint_mean=[x_u_i;x_u_j];
    joint_cov=blkdiag(P_u_i,P_u_j);
    
    [x_u_ij,P_u_ij]=linear_kf_update(joint_mean,joint_cov,A,b,Omega,R,z);
        
    x_joint(:,i)=x_u_ij;
    P_joint(:,:,i)=P_u_ij;
    
    if(marginal_calc(node_i)==0)
        x_u(:,node_i)=x_u_ij(1:2);
        P_u(:,:,node_i)=P_u_ij(1:2,1:2);
        marginal_calc(node_i)=1;
    end
    
    if(marginal_calc(node_j)==0)
        x_u(:,node_j)=x_u_ij(3:4);
        P_u(:,:,node_j)=P_u_ij(3:4,3:4);
        marginal_calc(node_j)=1;
    end
end

