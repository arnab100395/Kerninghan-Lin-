#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int nodes,edges,G=0,h=0,A[1000],v=0,k1=0, n1[1000][1000],w=0;
void d_value(int x[][nodes+1],int p1[],int p2[]);
void g_value(int a[][nodes+1],int d[],int p1[],int p2[]);
void swap1(int x, int y,int p1[], int p2[], int a[][nodes+1]);

int main(int argc, char *argv[])
{
    int v1, v2, i, j;
    for(i=0; i<100; i++)
    {
        for(j=0; j<100; j++)
        {
            n1[i][j]=0;
        }
    }
    for(i=0; i<100; i++)
    A[i]=0;
    FILE *infile;
    infile = fopen(argv[1],"r");
    if(infile==NULL)
    {
        printf("Not a valid filename");
        exit(0);
    }
    fscanf(infile,"%d",&nodes);
    w=nodes;
    if(nodes%2!=0)
    {
        printf("\nPlease enter even number of nodes\n");
        exit(0);
    }
    fscanf(infile,"%d",&edges);
    //printf("%d %d \n",nodes,edges);
    printf("\nInitial Solution: \n\n");
    printf("Part A  Part B\n");
    int edge[nodes+1][nodes+1],k, t=nodes/2,part1[t],part2[t];
    for(i=0;i<nodes+1;i++)
        {
        for(j=0;j<nodes+1;j++)
            edge[i][j]=0;
        }
    while(!feof(infile))
        {
            fscanf(infile,"%d",&v1);
            fscanf(infile,"%d",&v2);
            edge[v1][v2]=1;
        }
    for(i=0;i<nodes+1;i++)
        {
            for(j=0;j<nodes+1;j++)
                {
                    if (edge[i][j]==1)
                        edge[j][i]=1;
                    else
                        edge[i][j]=0;
                }
        }
    for( k=0;k<t;k++)
        {
            part1[k]=k+1;
            part2[k]=k+t+1;
        }
    for(i=0;i<t;i++)
        printf("%d       %d \n",part1[i],part2[i]);
    d_value(edge,part1,part2);
    return 0;
}

void d_value(int x[][nodes+1],int p1[],int p2[])
{
    int c[nodes+1],k,i,j;
    for(i=0;i<nodes+1;i++)
        {
            int e=0;
            for(j=0;j<nodes+1;j++)
                {
                    if(x[i][j]==1)
                    e++;
                }
            c[i]=e;
        }
    int t=nodes/2,l,d[nodes+1];
    //printf("\nD values \n");
    for(k=0;k<t;k++)
        {
            int cnt_c=0;
            for(l=0;l<t;l++)
                {
                    if(x[p1[k]][p2[l]]==1)
                    cnt_c++;
                }
            if (p1[k]==k+1)
                d[k+1]=2*cnt_c-c[k+1];
            else
                d[k+1]=0;
            //printf("%d ",d[k+1]);
        }
    //printf("\n");
    for(k=0;k<t;k++)
        {
            int cnt_c=0;
            for(l=0;l<t;l++)
                {
                    if(x[p2[k]][p1[l]]==1)
                    cnt_c++;
                }
            if(p2[k]==k+t+1)
                d[k+1+t]=2*cnt_c-c[k+1+t];
            else
                d[k+1+t]=0;
            //printf("%d ",d[k+1+t]);
        }
    g_value(x,d,p1,p2);
}

void g_value(int a[][nodes+1],int d[],int p1[],int p2[])
{
    int i,j,t=nodes/2,g[t+1][t+1],x=0,y=0,g_max=0,m=0,p;
    for (i=0;i<=t;i++)
        {
            for(j=0;j<=t;j++)
            g[i][j]=0;
        }
    for (i=1;i<=t;i++)
        {
            for(j=t+1;j<=nodes;j++)
                {
                    if (d[i]!=0 && d[j]!=0)
                        {
                            if (a[i][j]==1)
                                g[i][j-t]=d[i]+d[j]-2;
                            else
                                g[i][j-t]=d[i]+d[j];
                        }
                }
        }
    //printf("\n");
    /*for(i=0;i<=t;i++)
        {
            for(j=0;j<=t;j++)
                printf("%d ",g[i][j]);
            printf("\n");
        }*/
    int max=0;
    for (i=1;i<=t;i++)
        {
            for(j=t+1;j<=nodes;j++)
                {
                    if (g[i][j-t]>max)
                        {
                            max=g[i][j-t];
                            x=i; y=j;
                        }
                }
        }
    if (max==0)
        {
            for (i=1;i<=t;i++)
                {
                    for(j=t+1;j<=nodes;j++)
                        {
                        if (max==0 && max<abs(g[i][j-t]))
                                {
                                    max=abs(g[i][j-t]);
                                    x=i; y=j;
                                }
                        else if (max!=0 && abs(g[i][j-t])!=0 && max>abs(g[i][j-t]))
                        {
                            max=abs(g[i][j-t]);
                                    x=i; y=j;
                        }
                        }
                }max=-max;
         }
    G=G+max;
    A[h]=G;
    nodes=w;
    //printf("Nodes: %d",nodes);
    h++;
    //printf("G: %d \n",G);
    //printf("%d %d %d \n",x,y,max);
    //printf("g_max: ");
    if(G<=0)
        {
            //printf("Nodes: %d",nodes);
            if(h==1)
            {
                printf("\nAlready optimal, no swaps or further cuts");
                printf("\nFinal Solution: \n\n");
                printf("Part A  Part B\n");
                for (i=0;i<(nodes/2);i++)
                {
                   for(j=0;j<nodes;j++)
                        {
                            if (a[p1[i]][p2[j]]==1)
                            printf("(%d %d)\n",p1[i],p2[j]);
                        }
                }
                exit(0);
            }
            for(i=0;i<h-1;i++)
                {
                    if (A[i]>=g_max)
                        {
                            g_max=A[i];
                            m=i+1;
                        }
                }

            /*for(p=0;p<((nodes/2)-1);p++)
                {
                    printf("\n");
                    for(q=0;q<(nodes/2);q++)
                        printf("%d %d \n",n1[p][q],n1[p][q+(nodes/2)]);
                }
            printf("\n");*/
            printf("\nFinal Solution: \n\n");
            printf("Part A  Part B\n");
            //int q=0;
            for(p=0;p<(nodes/2);p++)
                printf("%d       %d \n",n1[m-1][p],n1[m-1][p+(nodes/2)]);
            printf("\nCut Set edges: \n");
            for (i=0;i<(nodes/2);i++)
                {
                   for(j=0;j<nodes;j++)
                        {
                            if (a[n1[m-1][i]][n1[m-1][j+(nodes/2)]]==1)
                            printf("(%d %d)\n",n1[m-1][i],n1[m-1][j+(nodes/2)]);
                        }
                }
            exit(0);
        }
    swap1(x,y,p1,p2,a);
}

void swap1(int x, int y,int p1[], int p2[], int a[][nodes+1])
{
    int  j, temp=0;
    temp=p1[x-1];
    p1[x-1]=p2[y-1-nodes/2];
    p2[y-1-nodes/2]=temp;

    for(j=0;j<(nodes/2);j++)
        {
            n1[k1][j]=p1[j];
            n1[k1][j+(nodes/2)]=p2[j];
        }
        //printf("\n");

    //for(j=0;j<nodes;j++)
      //  printf("%d ",n1[k1][j]);
    k1++;
    d_value(a,p1,p2);
}

