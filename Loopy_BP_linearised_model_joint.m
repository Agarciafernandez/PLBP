function [x_u,P_u,x_u_joint,P_u_joint]=Loopy_BP_linearised_model_joint(lin_parameters,edge_matrix,x_p,P_p,z_array,R,Nsteps)

%This function performs loopy belief propagation for the linearised model

%Author: Angel F. Garcia Fernandez

%Initialisation of messages
messages=struct;
messages_update=struct;

N=size(edge_matrix,2);


for k=1:Nsteps
    
    for i=1:N %We go through all the nodes, they transmit in parallel
        %Now we go through all neighbours
        list_neighbours_i=find(edge_matrix(i,:));
        [m_messages,n_messages]=size(messages);
        
        for j=1:length(list_neighbours_i)
            node_j=list_neighbours_i(j);
            list_neighbours_input_ij=list_neighbours_i([1:j-1,j+1:length(list_neighbours_i)]);
            
            %We first obtain the prior
            x_p_i=x_p(:,i);
            P_p_i=P_p(:,:,i);
            
            
            %We compute the KF updates
            
            x_u_i=x_p_i;
            P_u_i=P_p_i;
            
            for p=1:length(list_neighbours_input_ij)
                
                node_p=list_neighbours_input_ij(p); %node_p: neighbour of transmitting node
                
                if(and(node_p<=m_messages,i<=n_messages))  %If we haven't transmitted a message already, it is uniform so we don't do anything
                    alfa_p=messages(node_p,i).alfa;
                    if(~isempty(alfa_p))
                        
                        Hmess_p=messages(node_p,i).Hmess;
                        gamma_p=messages(node_p,i).gamma;
                        [x_u_i,P_u_i]=linear_kf_update(x_u_i,P_u_i,Hmess_p,0,0,gamma_p,alfa_p);
                    end
                    
                end
            end
            
            %We recover the linearisation parameters
            if(node_j<i)
                %Linearisation was done the other way
                A=lin_parameters(node_j,i).A;
                A_j=A(:,1:2);
                A_i=A(:,3:4);
                b=lin_parameters(node_j,i).b;
                Omega=lin_parameters(node_j,i).Omega;
                
            else
                A=lin_parameters(i,node_j).A;
                A_i=A(:,1:2);
                A_j=A(:,3:4);
                b=lin_parameters(i,node_j).b;
                Omega=lin_parameters(i,node_j).Omega;
            end
            
            
            alfa_sent=squeeze(z_array(i,node_j,:))-A_i*x_u_i-b;
            Hmess_sent=A_j;
            gamma_sent=Omega+R+A_i*P_u_i*A_i';
            
            messages_update(i,node_j).alfa=alfa_sent;
            messages_update(i,node_j).Hmess=Hmess_sent;
            messages_update(i,node_j).gamma=gamma_sent;
            
            
        end
        
        
    end
    messages=messages_update;
    

    
    
    
end

%Final calculation of marginal

[x_u_joint,P_u_joint,x_u,P_u]=jointMarginalCalculation(messages,lin_parameters,z_array,edge_matrix,x_p,P_p,R);



