clc;
clear;
Graph_list=graphofnet;
e=numnodes(Graph_list);
disp(e);
k1=ceil(e/2);
k2=floor((e-k1)*2/3);
k3=e-(k1+k2);
Adjacent_graph=adjacency(Graph_list);
cost=0;
add=0;
part_1=1:k1;  
part_2=[k1+1:k1+k2,zeros(1,k1-length(k1+1:k1+k2))];
part_3=[k1+k2+1:e,zeros(1,k1-length(k1+k2+1:e))];
matrix_initial=[part_1; part_2; part_3];
disp(matrix_initial);
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp('Initial partition');
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp((nonzeros(part_1))');
disp((nonzeros(part_2))');
disp((nonzeros(part_3))');
for i=1:50
    matrix_initial=[part_1; part_2; part_3];
    rand_select=randperm(3);
    m1=matrix_initial(rand_select(1),:);
    m2=matrix_initial(rand_select(2),:);
    [matrix_updated_1,matrix_updated_2,G]=twpp321(m1,m2,Adjacent_graph,k3);
    if isequal(m1,part_1)
        part_1=matrix_updated_1;
    elseif isequal(m1,part_2)
        part_2=matrix_updated_1;
    elseif isequal(m1,part_3)
        part_3=matrix_updated_1;
    end
    if isequal(m2,part_1)
        part_1=matrix_updated_2;
    elseif isequal(m2,part_2)
        part_2=matrix_updated_2;
    elseif isequal(m2,part_3)
        part_3=matrix_updated_2;
    end
    order=nonzeros([part_1 part_2 part_3 ]);
    Adjacent_graph=Adjacent_graph(order,order);
    part_2 = [part_2,zeros(1,k1-length(part_2))];  
    part_3 = [part_3,zeros(1,k1-length(part_3))]; 
    add=add+sum(sum(G));
    cost=[cost add];
    disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
    fprintf('\t\tPartition after Iteration Number %d : \n',i);
    disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
    disp((nonzeros(part_1))');
    disp((nonzeros(part_2))');
    disp((nonzeros(part_3))');
end
iter=0:i;
plot(iter,cost);
title('3-2-1 way weighted KL partitioning');
xlabel('Number of iterations');
ylabel('Cut Cost');
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp('Final Partition List :');
disp('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
disp((nonzeros(part_1))');
disp((nonzeros(part_2))');
disp((nonzeros(part_3))');
%Function to create the graph from the netlist
function Graph_list = graphofnet
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
%KL Function
function [final1,final2,g]= twpp321(matrix_1,matrix_2,A,k3)
n1=k3;
n2=k3;
n=n1+n2;
part1=matrix_1;
part2=matrix_2;
iteration=1;
done=0;
while ~done
    if iteration==1
    E_nc=[sum(A(1:n1,1:n1)) sum(A(n1+1:n,n1+1:n))];
    E_c=[sum(A(1:n1,n1+1:n)) sum(A(n1+1:n,1:n1))];
    D=E_c-E_nc;
    end
    gain=D(1:n1)'*ones(1,n2)+ ones(n1,1)*D(n1+1:n)-2*A(1:n1,n1+1:n);
    [temp,radd]=max(gain);
    [g(iteration),cadd]=max(temp);
    finalrowadd=radd(cadd); 
    finalcoladd=n1+cadd;
    swapnodeA(iteration)=part1(finalrowadd);
    swapnodeB(iteration)=part2(cadd);
    if size(A,1)==2
        done=1; % stop the next iteration, done
    else
        iteration=iteration+1; 
        part1=setdiff(part1,swapnodeA);
        part2=setdiff(part2,swapnodeB);
        idv=[setdiff([1:n],[finalrowadd finalcoladd]) [finalrowadd finalcoladd]]; % permute the indices of c
        A1 = A(idv,idv);
        D=D(idv(1:n-2));
        D(1:n1-1)=D(1:n1-1)+2*A1(n-1,1:n1-1)-2*A1(n,1:n1-1);
        D(n1:n-2)=D(n1:n-2)+2*A1(n,n1:n-2)-2*A1(n-1,n1:n-2);
        A= A1(1:n-2,1:n-2);
        n=n-2;
        n1=n1-1;
        n2=n2-1;
    end
end
% Evaluate the maximum Gain using the toeplitz function
[~,K]=max(sum(triu(toeplitz(g))));
%swap the nodes corresponding to the maximum gain
interm1=setdiff(matrix_1,swapnodeA(1:K));
final1=(nonzeros(union(interm1,swapnodeB(1:K))))';
interm2=setdiff(matrix_2,swapnodeB(1:K));
final2=(nonzeros(union(interm2,swapnodeA(1:K))))';
end

