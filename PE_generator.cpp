#include<bits/stdc++.h>
using namespace std;
struct atom{
    int x;
    int y;
    int z;
};
const double Navogadro=6.022e23;
int cell;
int space[505][505][505],self[505][505][505],mod[505][505][505];
atom pe[2080][5005],realpe[2080][5005];
atom trace[5100];
double dis(int n,int i,int j)
{
    double xx,yy,zz;
    xx=(pe[n][i].x-pe[n][j].x)*0.88912;
    yy=(pe[n][i].y-pe[n][j].y)*0.88912;
    zz=(pe[n][i].z-pe[n][j].z)*0.88912;
    double ans;
    ans=sqrt(pow(xx,2)+pow(yy,2)+pow(zz,2));
    return ans;
}
void selfpos()
{
    for(int i=0;i<2*cell+2;i++)
    {
        for(int j=0;j<2*cell+2;j++)
        {
            for(int k=0;k<2*cell+2;k++)
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
int natom=0,chain,chain_length[2080],length;
int main()
{
    int bonds=0,angles=0,dihedrals=0;
    long long res=0;
    double rho;
    freopen("chain_length.txt","r",stdin);
    string x;
    getline(cin,x);
    cin>>chain;
    srand(time(0));
    //fclose(stdin);
    //cin.clear();
    //freopen("muti_guass_distribution.txt","r",stdin);
    for(int i=0;i<chain;i++)
    {
        cin>>length;
        chain_length[i]=length;
        cout<<length;
        natom+=chain_length[i];
        if(chain_length[i]>=2) bonds+=chain_length[i]-1;
        if(chain_length[i]>=3) angles+=chain_length[i]-2;
        if(chain_length[i]>=4) dihedrals+=chain_length[i]-3;
    }
    //cell=ceil(pow(natom/1.5,1/3));
    //vector<vector<vector<int>>>space(cell*2,vector<vector<int>>(cell*2,vector<int>(cell*2,1)));
    cell=76;
    sort(chain_length,chain_length+chain,tmp);
    for(int i=0;i<2*cell+2;i++)//space.size()==2*cell
    {
        for(int j=0;j<2*cell+2;j++)
        {
            for(int k=0;k<2*cell+2;k++)
            {
                int kk=k%4;
                if(kk==0)
                {
                    if((i%4==0 && j%4==0 )||( j%4==2 && i%4==2))
                    space[i][j][k]=0;
                    else
                    space[i][j][k]=-1;
                }
                else if(kk==1)
                {
                    if(((i-1)%4==0 && (j-1)%4==0 )|| ((i-1)%4==2 && (j-1)%4==2))
                    space[i][j][k]=0;
                    else
                    space[i][j][k]=-1;
                }
                else if(kk==2)
                {
                    if((i%4==0 && j%4==2 )|| (i%4==2 && j%4==0))
                    space[i][j][k]=0;
                    else
                    space[i][j][k]=-1;
                }
                else
                {
                    if((i%4==3 && j%4==1 )||( i%4==1 && j%4==3))
                    space[i][j][k]=0;
                    else
                    space[i][j][k]=-1;
                }
                
            }
        }
    }
    for(int i=0;i<2*cell+2;i++)
    {
        for(int j=0;j<2*cell+2;j++)
        {
            for(int k=0;k<2*cell+2;k++)
            {
                mod[i][j][k]=space[i][j][k];
            }
        }
    }

    int xpos,ypos,zpos,xx,yy,zz,flag=1;
    int rx[2500],ry[2500],rz[2500];
    int step_ud,step_lr;//up&down,left&right;
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
                    xpos=(rand() % (2*cell));
                    ypos=(rand() % (2*cell));
                    zpos=(rand() % (2*cell));
                }while(space[xpos][ypos][zpos]==1);
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
                int t=0;
                flag=1;
                do{
                    t++;
                    if(t==100)   //四个方向都受阻，此链重开；
                    {
                        flag=0;
                        goto LABEL1;
                    }
                    xx=xpos;
                    yy=ypos;
                    zz=zpos;//xx,yy,zz复制父节点的位置
                    step_ud=pow(-1,(rand() % 2)+1);//-1为下，1为上
                    step_lr=pow(-1,(rand() % 2)+1);//-1为左，1为右
                    rz[i]=rz[i-1]+step_ud;//rz real z
                    zz+=step_ud;
                    if(zz>=2*cell) zz-=2*cell;
                    if(zz<=0)zz+=2*cell;
                    if(zpos%2==0)//偶数层
                    {
                        if(step_ud==-1)//偶数层，向下走
                        {
                            xx+=step_lr;
                            yy-=step_lr;
                            rx[i]=rx[i-1]+step_lr;
                            ry[i]=ry[i-1]-step_lr;
                        }
                        else
                        {
                            xx+=step_lr;
                            yy+=step_lr;
                            rx[i]=rx[i-1]+step_lr;
                            ry[i]=ry[i-1]+step_lr;
                        }
                    }
                    else//奇数层
                    {
                        if(step_lr==-1)//奇数层，向下走
                        {
                            xx+=step_lr;
                            yy+=step_lr;
                            rx[i]=rx[i-1]+step_lr;
                            ry[i]=ry[i-1]+step_lr;
                        }
                        else
                        {
                            xx+=step_lr;
                            yy-=step_lr;
                            rx[i]=rx[i-1]+step_lr;
                            ry[i]=ry[i-1]-step_lr;
                        }
                    }
                    
                    if(xx>=2*cell)xx-=2*cell;
                    if(xx<=0)xx+=2*cell;
                    if(yy>=2*cell)yy-=2*cell;
                    if(yy<=0)yy+=2*cell;//zz,xx,yy在复制的父节点基础上走动一格。
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
                        if(i!=chain_length[n])
                        goto LABLE2;
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
    rho=natom/Navogadro*14.02/pow(2*cell*0.88913*1e-8,3);
    //-----------生成链的周期性边界位置---------
    freopen("wrapped_coordinate.dat","w",stdout);
    cout<<"#MODEL FOR PE WITH A DISTRIBUTION OF SCHULTZ-ZIMM  "<<res<<endl;
    cout<<"#SYSTEM DENSITY:"<<rho<<"g/cm3"<<endl;
    cout<<natom<<"  "<<"atoms"<<endl;
    cout<<bonds<<" "<<"bonds"<<endl;//加和 每条链的长度减一
    cout<<angles<<" "<<"angles"<<endl;
    cout<<dihedrals<<" "<<"dihedrals"<<endl;
    cout<<endl<<"1     atom types"<<endl<<"1     bond types"<<endl<<"1     angle types"<<endl<<"1     dihedral types"<<endl;
    cout<<"0"<<"    "<<2*cell*0.88913<<" xlo xhi"<<endl;
    cout<<"0"<<"    "<<2*cell*0.88913<<" ylo yhi"<<endl;
    cout<<"0"<<"    "<<2*cell*0.88913<<" zlo zhi"<<endl;
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
            cout<<num<<"    "<<n+1<<"    "<<"1"<<"   "<<pe[n][i].x*0.88912<<"    "<<pe[n][i].y*0.88912<<"    "<<pe[n][i].z*0.88912<<endl;
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
    freopen("unwrapped_coordinate.dat","w",stdout);
    cout<<"#MODEL FOR PE WITH DIFFERENT DP  "<<res<<endl;
    cout<<natom<<"  "<<"atoms"<<endl;
    cout<<bonds<<" "<<"bonds"<<endl;//加和 每条链的长度减一
    cout<<angles<<" "<<"angles"<<endl;
    cout<<dihedrals<<" "<<"dihedrals"<<endl;
    cout<<endl<<"1     atom types"<<endl<<"1     bond types"<<endl<<"1     angle types"<<endl<<"1     dihedral types"<<endl;
    cout<<"0"<<"    "<<2*cell*0.88912<<" xlo xhi"<<endl;
    cout<<"0"<<"    "<<2*cell*0.88912<<" ylo yhi"<<endl;
    cout<<"0"<<"    "<<2*cell*0.88912<<" zlo zhi"<<endl;
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
            cout<<num<<"    "<<n+1<<"    "<<"1"<<"   "<<realpe[n][i].x*0.88912<<"    "<<realpe[n][i].y*0.88912<<"    "<<realpe[n][i].z*0.88912<<endl;
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