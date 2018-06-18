% main program
clc;
clear;
Graph_list=graph_net;
nodes=mod(numnodes(Graph_list),5);
if(nodes~=0)
    return;
end
nodes=(numnodes(Graph_list))/5;
adjacent_graph=adjacency(Graph_list);  
cost=0;
add=0;
part_1=0*nodes+1:1*nodes;
part_2=1*nodes+1:2*nodes;
part_3=2*nodes+1:3*nodes;
part_4=3*nodes+1:4*nodes;
part_5=4*nodes+1:5*nodes;
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp('Initial partition');
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp(part_1);
disp(part_2);
disp(part_3);
disp(part_4);
disp(part_5);
for i=1:100
    matrix_initial=[part_1; part_2; part_3; part_4; part_5];
    random_select=randperm(5);
    m1=matrix_initial(random_select(1),:);
    m2=matrix_initial(random_select(2),:);
    [matrix_updated_1,matrix_updated_2,G]=twpp(m1,m2,adjacent_graph,nodes);
    if (m1==part_1)
        part_1=matrix_updated_1;
    elseif (m1==part_2)
        part_2=matrix_updated_1;
    elseif (m1==part_3)
        part_3=matrix_updated_1;
    elseif (m1==part_4)
        part_4=matrix_updated_1;
    else
        part_5=matrix_updated_1;
    end
    if (m2==part_1)
        part_1=matrix_updated_2;
    elseif (m2==part_2)
        part_2=matrix_updated_2;
    elseif (m2==part_3)
        part_3=matrix_updated_2;
    elseif (m2==part_4)
        part_4=matrix_updated_2;
    else
        part_5=matrix_updated_2;
    end
    order=[part_1 part_2 part_3 part_4 part_5];
    adjacent_graph=adjacent_graph(order,order);
    add=add+sum(sum(G));
    cost=[cost add];
    disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
    fprintf('\t\tPartition after Iteration Number %d : \n',i);
    disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
    disp(part_1);
    disp(part_2);
    disp(part_3);
    disp(part_4);
    disp(part_5);
end
iter=0:i;
plot(iter,cost);
title('5 way even KL partitioning');
xlabel('Number of iterations');
ylabel('Cut Cost');
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp('Final Partition List :');
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp(part_1);
disp(part_2);
disp(part_3);
disp(part_4);
disp(part_5);

% function converting a netlist to a graph
function Graph_list = graph_net
clc;
clear;
row_size_max =1000;
col_size_max=100;
m=1;
file1= fopen('net.txt','r'); 
while ~feof(file1)
    line = textscan(file1,'%s%s',row_size_max,'Delimiter',':');
end
new_array=zeros(row_size_max,col_size_max);
for i=1:size(line{2})
    line1=sscanf(line{2}{i},'%u');
    for j=1:size(line1)
        new_array(i,j)=line1(j);
    end
end
size_max=max(new_array(:));
array_final=zeros((size_max-1)*i,2);
for i=1:size(line{2})
    for j=1:size_max-1
        if (new_array(i,j)~=0 && new_array(i,j+1)~=0)
            array_final(m,1)=new_array(i,j);
            array_final(m,2)=new_array(i,j+1);
            m=m+1;
        end
    end
end
unique_array = unique ([min(array_final(:,[1,2]),[],2) max(array_final(:,[1,2]),[],2)],'rows');
unique_array(1,:)=[];
fclose(file1);
row_number = size(unique_array,1); 
edge_weights=(ones(1,row_number))';
EdgeTable=table(unique_array,edge_weights,'VariableNames',{'EndNodes' 'Weight'});
Graph_list = graph(EdgeTable);
end
%d A function for a 2 way Kerninghan Lin Partition
function [final1,final2,G]= twpp(partition1,partition2,Adjacent_graph,nodes)
n1=nodes;
n2=nodes;
n=n1+n2;
part1=partition1;
part2=partition2;
iteration=1;
done=0;
E_not_connected=[sum(Adjacent_graph(1:n1,1:n1)) sum(Adjacent_graph(n1+1:n,n1+1:n))];
E_connected=[sum(Adjacent_graph(1:n1,n1+1:n)) sum(Adjacent_graph(n1+1:n,1:n1))];
D=E_connected-E_not_connected;
while ~done
    gain=D(1:n1)'*ones(1,n2)+ ones(n1,1)*D(n1+1:n)-2*Adjacent_graph(1:n1,n1+1:n);
    [temp,radd]=max(gain);
    [g(iteration),cadd]=max(temp);
    finalrowadd=radd(cadd); 
    swapnodeA(iteration)=part1(finalrowadd);
    swapnodeB(iteration)=part2(cadd);
    if size(Adjacent_graph,1)==2
        done=1;
    else
    iteration=iteration+1; 
    part1=setdiff(part1,swapnodeA);
    part2=setdiff(part2,swapnodeB);
    D(1:n1-1)=D(1:n1-1)+2*Adjacent_graph(n-1,1:n1-1)-2*Adjacent_graph(n,1:n1-1);
    D(n1:n-2)=D(n1:n-2)+2*Adjacent_graph(n,n1:n-2)-2*Adjacent_graph(n-1,n1:n-2);
    Adjacent_graph= Adjacent_graph(1:n-2,1:n-2);
    n=n-2;
    n1=n1-1;
    n2=n2-1;
    end
end
G=g;
% Gain is evaluated using the Toeplitz function
[~,K]=max(sum(triu(toeplitz(g))));
% Nodes are now swapped 
interm1=setdiff(partition1,swapnodeA(1:K));
final1=union(interm1,swapnodeB(1:K));
interm2=setdiff(partition2,swapnodeB(1:K));
final2=union(interm2,swapnodeA(1:K));
end