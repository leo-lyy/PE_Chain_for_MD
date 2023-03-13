#include<bits/stdc++.h>
using namespace std;
struct atom{
    int x;
    int y;
    int z;
};
const double Navogadro=6.022e23;
int cell,chain,length,chain_length[2080];
long int natom=0;
int space[505][505][505],self[505][505][505],mod[505][505][505];
atom pe[2080][5005],realpe[2080][5005];
atom trace[5100];
void selfpos()
{
    for(int i=0;i<cell+2;i++)
    {
        for(int j=0;j<cell+2;j++)
        {
            for(int k=0;k<cell+2;k++)
            {
                self[i][j][k]=mod[i][j][k];
            }
        }
    }

}
bool tmp(int a,int b)
{
    return a>b;
}
int main()
{
    int bonds=0,angles=0,dihedrals=0;
    long long res=0;
    double rho;
    freopen("chain_length.txt","r",stdin);
    string x;
    getline(cin,x);
    cin>>chain;
    cell=500;
    srand(time(0));
    int lennum;
    for(int i=0;i<chain;i++)
    {
        cin>>length;
        cin>>lennum;
        for (int k=0;k<lennum;k++)
        {
            chain_length[i+k]=length;
            //cout<<length;
            natom+=chain_length[i+k];
            if(chain_length[i+k]>=2) bonds+=chain_length[i+k]-1;
            if(chain_length[i+k]>=3) angles+=chain_length[i+k]-2;
            if(chain_length[i+k]>=4) dihedrals+=chain_length[i+k]-3;
        }
        i+=lennum-1;

    }
    fclose(stdin);
    freopen("CON", "r", stdin);
    char ch;
    while(1)
    {
        cout<<"Please input a box size: "<<endl;
        cin>>cell;
        rho=natom/Navogadro*14.02/pow(cell*1.082*1e-8,3);
        cout<<"the density of the box is "<<rho<<"g/cm3"<<endl;
        cout<<"continue?(Y/n)";
        cin>>ch;
        if(ch=='Y'||ch=='y')
        {
            break;
        }
        else continue;
            
    }
    int i,j,k;
    for(k=0;k<cell;k++)
    {
        for(i=0;i<cell;i++)
        {
            for(j=0;j<cell;j++)
            {
                if(k%2==0)
                {
                   space[i][j][k]=(i+j)%2-1; //非fcc点 -1
                   mod[i][j][k]=space[i][j][k];
                }
                else
                {
                    space[i][j][k]=(i+j+1)%2-1; //fcc晶格点 0
                    mod[i][j][k]=space[i][j][k];
                }
            }
        }
    }
    //-------------------------------

   
    int xpos,ypos,zpos,xx,yy,zz,flag=1;
    int rx[2500],ry[2500],rz[2500];
    for(int n=0;n<chain;)
    {
        LABEL1: selfpos();
        res++;
        cout<<"Generating chain    "<<n+1<<endl;
        for(int i=0;i<chain_length[n];i++)
        {
            if(i==0)
            {
                do{
                    xpos=(rand() % (cell-2))+1;
                    ypos=(rand() % (cell-2))+1;
                    zpos=(rand() % (cell-2))+1;
                }while(space[xpos][ypos][zpos]!=0);
                //space[xpos][ypos][zpos]=1;
                trace[i].x=xpos;
                trace[i].y=ypos;
                trace[i].z=zpos;
                self[xpos][ypos][zpos]=1;
                rx[0]=xpos;
                ry[0]=ypos;
                rz[0]=zpos;
                realpe[n][i].x=rx[0];
                realpe[n][i].y=ry[0];
                realpe[n][i].z=rz[0];
            }
            else 
            {
                int t=0,step_ud,step_locate;
                flag=1;
                do{
                    t++;
                    if(t==1200)   //12个方向都受阻，此链重开；
                    {
                        flag=0;
                        goto LABEL1;
                    }
                    xx=xpos;
                    yy=ypos;
                    zz=zpos;//xx,yy,zz复制父节点的位置
                    step_ud=(rand()%3)-1;//-1为向下，1为向上，0为当前层

                    step_locate=rand()%4; //0,1,2,3

                    rz[i]=rz[i-1]+step_ud;//rz real z
                    zz+=step_ud;
                    if(zz>=cell) zz-=cell;
                    if(zz<=0) zz+=cell;

                    if(step_ud%2!=0)//上下走
                    {
                        if(step_locate==0)
                        {
                            yy+=1;
                            rx[i]=rx[i-1];
                            ry[i]=ry[i-1]+1;
                        }
                        if(step_locate==1)
                        {
                            xx-=1;
                            rx[i]=rx[i-1]-1;
                            ry[i]=ry[i-1];
                        }
                        if(step_locate==2)
                        {
                            yy-=1;
                            rx[i]=rx[i-1];
                            ry[i]=ry[i-1]-1;
                        }
                        if(step_locate==3)
                        {
                            xx+=1;
                            rx[i]=rx[i-1]+1;
                            ry[i]=ry[i-1];
                        }
                    }
                    else//当前层
                    {
                        if(step_locate==0)
                        {
                            xx+=1;
                            yy+=1;
                            rx[i]=rx[i-1]+1;
                            ry[i]=ry[i-1]+1;
                        }
                        if(step_locate==1)
                        {
                            xx-=1;
                            yy+=1;
                            rx[i]=rx[i-1]-1;
                            ry[i]=ry[i-1]+1;
                        }
                        if(step_locate==2)
                        {
                            xx-=1;
                            yy-=1;
                            rx[i]=rx[i-1]-1;
                            ry[i]=ry[i-1]-1;
                        }
                        if(step_locate==3)
                        {
                            xx+=1;
                            yy-=1;
                            rx[i]=rx[i-1]+1;
                            ry[i]=ry[i-1]-1;
                        }
                    }
                    
                    if(xx>=cell)xx-=cell;
                    if(xx<=0)xx+=cell;
                    if(yy>=cell)yy-=cell;
                    if(yy<=0)yy+=cell;//zz,xx,yy在复制的父节点基础上走动一格。
                    if(space[xx][yy][zz]==0&&self[xx][yy][zz]==0)//走动到的新节点没有被空间其他链或自身链占据
                    {
                        self[xx][yy][zz]=1; //自身链新增的空间坐标置为1
                        xpos=xx;
                        ypos=yy;
                        zpos=zz;//父节点更新
                        trace[i].x=xpos;
                        trace[i].y=ypos;
                        trace[i].z=zpos;//记录自身链新增原子的周期性坐标位置
                        // realpe[n][i].x=rx;
                        // realpe[n][i].y=ry;
                        // realpe[n][i].z=rz;//记录自身链新增原子的真实位置
                        if(i!=chain_length[n]) goto LABLE2;
                        break;
                    }
                    else continue;
                }while(1);
                if(flag==1)
                {
                    
                    //self[xpos][ypos][zpos]=1;
                }
                else 
                {
                    break;
                }
                
            } 
            LABLE2 :continue;
        }
        if(flag==1)
        {
            selfpos();
            for(int i=0;i<chain_length[n];i++)
            {
                pe[n][i].x=trace[i].x;
                pe[n][i].y=trace[i].y;
                pe[n][i].z=trace[i].z;
                realpe[n][i].x=rx[i];
                realpe[n][i].y=ry[i];
                realpe[n][i].z=rz[i];
                space[trace[i].x][trace[i].y][trace[i].z]=1;
            }
            n++;
        }

    }

int num=0,len=0;
    rho=natom/Navogadro*14.02/pow(cell*1.082*1e-8,3);
    //-----------生成链的周期性边界位置---------
    freopen("fcc_wrapped_coordinate.dat","w",stdout);
    cout<<"#MODEL FOR PE WITH A DISTRIBUTION OF SCHULTZ-ZIMM  "<<res<<endl;
    cout<<"#SYSTEM DENSITY:"<<rho<<"g/cm3"<<endl;
    cout<<natom<<"  "<<"atoms"<<endl;
    cout<<bonds<<" "<<"bonds"<<endl;//加和 每条链的长度减一
    cout<<angles<<" "<<"angles"<<endl;
    cout<<dihedrals<<" "<<"dihedrals"<<endl;
    cout<<endl<<"1     atom types"<<endl<<"1     bond types"<<endl<<"1     angle types"<<endl<<"1     dihedral types"<<endl;
    cout<<"0"<<"    "<<cell*1.082<<" xlo xhi"<<endl;
    cout<<"0"<<"    "<<cell*1.082<<" ylo yhi"<<endl;
    cout<<"0"<<"    "<<cell*1.082<<" zlo zhi"<<endl;
    cout<<endl<<"Masses"<<endl;
    cout<<endl<<"1          14.02"<<endl;
    cout<<endl<<"Atoms"<<endl;
    cout<<endl;

    //cout<<chain_length<<endl;
    for(int n=0;n<chain;n++)
    {
        for(int i=0;i<chain_length[n];i++)
        {
            num++;
            cout<<num<<"    "<<n+1<<"    "<<"1"<<"   "<<pe[n][i].x*1.082<<"    "<<pe[n][i].y*1.082<<"    "<<pe[n][i].z*1.082<<endl;
           // cout<<num<<"    "<<n+1<<"    "<<"1"<<"   "<<pe[n][i].x<<"    "<<pe[n][i].y<<"    "<<pe[n][i].z<<endl;
        
        }
        //cout<<chain_length[n]<<endl;
    }
    num=0;
    cout<<endl<<"Bonds"<<endl<<endl;
    for(int n=0;n<chain;n++)
    {
        if(n>0)len+=chain_length[n-1];
        for(int i=0;i<chain_length[n]-1;i++)
        {
            num++;
            if(chain_length[n]<2)break;
            cout<<num<<"    "<<"1"<<"   "<<i+1+len<<"    "<<i+2+len<<endl;///////
        }
    }
    cout<<endl<<"Angles"<<endl<<endl;
    num=0;
    len=0;
    for(int n=0;n<chain;n++)
    {
        if(n>0)len+=chain_length[n-1];
        for(int i=0;i<chain_length[n]-2;i++)
        {
            num++;
            if(chain_length[n]<3)break;
            cout<<num<<"    "<<"1"<<"   "<<i+1+len<<"    "<<i+2+len<<"    "<<i+3+len<<endl;
        }
    }
    cout<<endl<<"Dihedrals"<<endl<<endl;
    num=0;
    len=0;
    for(int n=0;n<chain;n++)
    {
        if(n>0)len+=chain_length[n-1];
        for(int i=0;i<chain_length[n]-3;i++)
        {
            num++;
            if(chain_length[n]<3)break;
            cout<<num<<"    "<<"1"<<"   "<<i+1+len<<"    "<<i+2+len<<"    "<<i+3+len<<"    "<<i+4+len<<endl;
        }
    }
    //-----------生成链的实际位置---------
    freopen("fcc_unwrapped_coordinate.dat","w",stdout);
    cout<<"#MODEL FOR PE WITH DIFFERENT DP  "<<res<<endl;
    cout<<natom<<"  "<<"atoms"<<endl;
    cout<<bonds<<" "<<"bonds"<<endl;//加和 每条链的长度减一
    cout<<angles<<" "<<"angles"<<endl;
    cout<<dihedrals<<" "<<"dihedrals"<<endl;
    cout<<endl<<"1     atom types"<<endl<<"1     bond types"<<endl<<"1     angle types"<<endl<<"1     dihedral types"<<endl;
    cout<<"0"<<"    "<<cell*1.082<<" xlo xhi"<<endl;
    cout<<"0"<<"    "<<cell*1.082<<" ylo yhi"<<endl;
    cout<<"0"<<"    "<<cell*1.082<<" zlo zhi"<<endl;
    cout<<endl<<"Masses"<<endl;
    cout<<endl<<"1          14.02"<<endl;
    cout<<endl<<"Atoms"<<endl;
    cout<<endl;

    //cout<<chain_length<<endl;
    num=0;
    for(int n=0;n<chain;n++)
    {
        for(int i=0;i<chain_length[n];i++)
        {
            num++;
            cout<<num<<"    "<<n+1<<"    "<<"1"<<"   "<<realpe[n][i].x*1.082<<"    "<<realpe[n][i].y*1.082<<"    "<<realpe[n][i].z*1.082<<endl;
        //cout<<num<<"    "<<n+1<<"    "<<"1"<<"   "<<realpe[n][i].x<<"    "<<realpe[n][i].y<<"    "<<realpe[n][i].z<<endl;
        
        }
        //cout<<chain_length[n]<<endl;
    }
    num=0;
    len=0;
    cout<<endl<<"Bonds"<<endl<<endl;
    for(int n=0;n<chain;n++)
    {
        if(n>0)len+=chain_length[n-1];
        for(int i=0;i<chain_length[n]-1;i++)
        {
            num++;
            cout<<num<<"    "<<"1"<<"   "<<i+1+len<<"    "<<i+2+len<<endl;
        }
    }
    cout<<endl<<"Angles"<<endl<<endl;
    num=0;
    len=0;
    for(int n=0;n<chain;n++)
    {
        if(n>0)len+=chain_length[n-1];
        for(int i=0;i<chain_length[n]-2;i++)
        {
            num++;
            cout<<num<<"    "<<"1"<<"   "<<i+1+len<<"    "<<i+2+len<<"    "<<i+3+len<<endl;
        }
    }
    cout<<endl<<"Dihedrals"<<endl<<endl;
    num=0;
    len=0;
    for(int n=0;n<chain;n++)
    {
        if(n>0)len+=chain_length[n-1];
        for(int i=0;i<chain_length[n]-3;i++)
        {
            num++;
            cout<<num<<"    "<<"1"<<"   "<<i+1+len<<"    "<<i+2+len<<"    "<<i+3+len<<"    "<<i+4+len<<endl;
        }
    }
    /*
    for(int k=0;k<2*cell-1;k++)
    {
        for(int j=0;j<2*cell;j++)
        {
            for(int i=0;i<2*cell;i++)
            {
                if(space[i][j][k]==0)
                cout<<"0"<<"    ";
                else
                cout<<"*"<<"    ";
            }
            cout<<endl;
        }
        cout<<"---------------------------------"<<endl;
    }
    */
   return 0;
}
