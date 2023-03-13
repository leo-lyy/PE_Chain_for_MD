#include<bits/stdc++.h>
using namespace std;
struct atom{
    int x;
    int y;
    int z;
};
struct mol{
    int bonelen;
    int branchnum;
    int branchpos[105]={0};
    int branchlen[105]={0};
};
const double Navogadro=6.022e23;
atom pe[2080][5005],realpe[2080][5005];  //realpe 实际位置
atom trace[5100];   //随机生成中的链行走轨迹
mol comb[2090];
int chain,natom,cell;
int space[505][505][505],self[505][505][505],mod[505][505][505];
int rx[2500],ry[2500],rz[2500];


void selfpos() //初始化当前生成链
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


int generateatom(int i,int f,int *px,int *py,int *pz) //生成第i个原子，i原子的父原子为f，后三个参数为此原子父节点坐标
{
    int t=0,step_ud,step_locate,xx,yy,zz;
    int flag=1;
    do{
        t++;
        if(t==1200)  //12个方向都受阻，此链重开；
        {
            flag=0;
            return 0;
        }
        xx= *px;
        yy= *py;
        zz= *pz;//xx,yy,zz复制父节点的位置
        step_ud=(rand()%3)-1;//-1为向下，1为向上，0为当前层
        step_locate=rand()%4; //0,1,2,3
        rz[i]=rz[f]+step_ud;//rz real z
        zz+=step_ud;
        if(zz>=cell) zz-=cell;
        if(zz<=0) zz+=cell;

        if(step_ud%2!=0)//上下走
        {
            if(step_locate==0)
            {
                yy+=1;
                rx[i]=rx[f];          //rx ry rz 实际原子坐标位置
                ry[i]=ry[f]+1;
            }
            if(step_locate==1)
            {
                xx-=1;
                rx[i]=rx[f]-1;
                ry[i]=ry[f];
            }
            if(step_locate==2)
            {
                yy-=1;
                rx[i]=rx[f];
                ry[i]=ry[f]-1;
            }
            if(step_locate==3)
            {
                xx+=1;
                rx[i]=rx[f]+1;
                ry[i]=ry[f];
            }
        }
        else//当前层
        {
            if(step_locate==0)
            {
                xx+=1;
                yy+=1;
                rx[i]=rx[f]+1;
                ry[i]=ry[f]+1;
            }
            if(step_locate==1)
            {
                xx-=1;
                yy+=1;
                rx[i]=rx[f]-1;
                ry[i]=ry[f]+1;
            }
            if(step_locate==2)
            {
                xx-=1;
                yy-=1;
                rx[i]=rx[f]-1;
                ry[i]=ry[f]-1;
            }
            if(step_locate==3)
            {
                xx+=1;
                yy-=1;
                rx[i]=rx[f]+1;
                ry[i]=ry[f]-1;
            }
        }
        
        if(xx>=cell)xx-=cell;
        if(xx<=0)xx+=cell;
        if(yy>=cell)yy-=cell;
        if(yy<=0)yy+=cell;//zz,xx,yy在复制的父节点基础上走动一格。
        if(space[xx][yy][zz]==0&&self[xx][yy][zz]==0)//走动到的新节点没有被空间其他链或自身链占据
        {
            //cout<<i<<": "<<xx<<" "<<yy<<" "<<zz<<endl;
            cout<<"#"<<i<<" "<<f<<": "<<rx[i]<<" "<<ry[i]<<" "<<rz[i]<<endl;
            self[xx][yy][zz]=1; //自身链新增的空间坐标置为1
            *px=xx;
            *py=yy;
            *pz=zz;//父节点更新
            trace[i].x=xx;
            trace[i].y=yy;
            trace[i].z=zz;//记录自身链新增原子的周期性坐标位置
            return 1;

            // if(i!=chain_length[n])
            // {
            //     return 1;
            // }
            // break;
        }
    }while(1);
}
void writechain(int n,int bac)
{
    for(int i=0;i<=comb[n].bonelen+bac;i++)
    {
        pe[n][i].x=trace[i].x;
        pe[n][i].y=trace[i].y;
        pe[n][i].z=trace[i].z;
        realpe[n][i].x=rx[i];
        realpe[n][i].y=ry[i];
        realpe[n][i].z=rz[i];
        space[trace[i].x][trace[i].y][trace[i].z]=1;
        selfpos();
    }

}

int main()
{
    int bonds=0,angles=0,dihedrals=0,sideatom=0;
    long long res=0;
    double rho;
    cell=100;
    natom=0;
    freopen("chain_length.txt","r",stdin);
    string x;
    getline(cin,x);
    cin>>chain;
    srand(time(0));
    for(int i=0;i<chain;i++)
    {
        cin>>comb[i].bonelen>>comb[i].branchnum;
        natom+=comb[i].bonelen;
        sideatom=0;
        for (int a=1;a<=comb[i].branchnum;a++)
        { 
            cin>>comb[i].branchpos[a]>>comb[i].branchlen[a];
            sideatom+=comb[i].branchlen[a];
            natom+=comb[i].branchlen[a];
        }
        //cout<<length;
        bonds+=(comb[i].bonelen+sideatom)-1;
        angles+=(comb[i].bonelen-2)+(sideatom-comb[i].branchnum)+2*comb[i].branchnum;
        dihedrals+=(comb[i].bonelen-3)+(sideatom-2*comb[i].branchnum)+comb[i].branchnum*4;
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

    int i,j,k;  //初始化fcc晶格
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
    int *px,*py,*pz;
    px= &xpos;
    py= &ypos;
    pz= &zpos;

    for(int n=0;n<chain;n++)
    {
        LABEL1: selfpos();
        res++;
        cout<<"Generating chain    "<<n+1<<endl;
        int b=1,bac=0,f;  //b支链序号，bac:branch atom count 
        for(int i=0;i<comb[n].bonelen;i++)
        {
            if(i==0)//链的第一个原子
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
            else if (i+1==comb[n].branchpos[b]) //原子为支化点
            {
                f=i-1;
                int gen=generateatom(i,f,px,py,pz); // 生成支化点原子
                if (gen==0)
                {
                    cout<<"chain regenerating 1..."<<endl;
                    goto LABEL1;
                }
                
                int bx=xpos,by=ypos,bz=zpos; //记录支化点的坐标
                int f=i;
                for (int k=1;k<=comb[n].branchlen[b];k++) //生成支链上的原子
                {
                    bac++;
                    int gen=generateatom(comb[n].bonelen+bac,f,px,py,pz);
                    if (gen==0)
                    {
                        cout<<"chain regenerating 2..."<<endl;
                        goto LABEL1;
                    }
                    f=comb[n].bonelen+bac;
                }
                b++;
                xpos=bx;
                ypos=by;
                zpos=bz;
            }
            else                                //生成主链上的原子
            {
                f=i-1;
                int gen=generateatom(i,f,px,py,pz);
                if (gen==0)
                {
                    cout<<"chain regenerating 3..."<<endl;
                    goto LABEL1;
                }
                else
                {
                    writechain(n,bac);
                }

            }
        }
        //当前链生成完成，将链数据写入map
        cout<<"------------------"<<endl;

    }
    // system("pause");
    freopen("fcc_unwrapped_coordinate.dat","w",stdout);
    cout<<"#MODEL FOR PE WITH DIFFERENT DP  "<<res<<endl;
    cout<<natom<<"  "<<"atoms"<<endl;
    cout<<natom-chain<<" "<<"bonds"<<endl;//加和 每条链的长度减一
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
    int bac=0,b=1;
    float sc=1.082;
    long int count=0,mol=0;
    for(int n=0;n<chain;n++)
    {
        bac=0,b=1;
        for(int i=0;i<comb[n].bonelen;i++)
        {
            if(i+1==comb[n].branchpos[b])
            {
                printf("%d\t%d\t1\t%.4f\t%.4f\t%.4f\n",count+i+1,n+1,sc*realpe[n][i].x,sc*realpe[n][i].y,sc*realpe[n][i].z);
                for (int j=1;j<=comb[n].branchlen[b];j++)
                {
                    bac++;
                    printf("%d\t%d\t1\t%.4f\t%.4f\t%.4f\n",count+comb[n].bonelen+bac,n+1,sc*realpe[n][comb[n].bonelen+bac].x,sc*realpe[n][comb[n].bonelen+bac].y,sc*realpe[n][comb[n].bonelen+bac].z);

                }
                b++;
            }
            else
            {
                printf("%d\t%d\t1\t%.4f\t%.4f\t%.4f\n",count+i+1,n+1,1.000*sc*realpe[n][i].x,1.000*sc*realpe[n][i].y,1.000*sc*realpe[n][i].z);
            }
        }
        count+=comb[n].bonelen+bac;
    }
    count=0;

    cout<<endl<<"Bonds"<<endl<<endl;
    for(int n=0;n<chain;n++)
    {
        bac=0,b=1;
        for(int i=0;i<comb[n].bonelen-1;i++)
        {
            if(i+1==comb[n].branchpos[b])
            {
                cout<<count+comb[n].bonelen+bac<<"  "<<"1"<<"   "<<count+i+1+mol<<"   "<<count+comb[n].bonelen+bac+1+mol<<endl;
                for(int j=0;j<comb[n].branchlen[b]-1;j++)
                {
                    cout<<count+comb[n].bonelen+bac+1<<"    "<<" "<<"1"<<"   "<<count+comb[n].bonelen+bac+2+mol<<"   "<<count+comb[n].bonelen+bac+1+mol<<endl;
                    bac++;
                }
                b++;
                bac++;
                cout<<count+i+1<<"  "<<"1"<<"   "<<count+i+1+mol<<"   "<<count+i+2+mol<<endl;
            }
            else
            {
                cout<<count+i+1<<"  "<<"1"<<"   "<<count+i+1+mol<<"   "<<count+i+2+mol<<endl;
            }
        }
        count+=comb[n].bonelen+bac-1;
        mol++;
    }
    count=0;

    cout<<endl<<"Angles"<<endl<<endl;
    for(int n=0;n<chain;n++)
    {
        bac=0,b=1;
        for(int i=1;i<comb[n].bonelen-1;i++)
        {
            if(i+1==comb[n].branchpos[b])
            {
                cout<<count+i<<"  "<<"1"<<"   "<<count+i<<"   "<<count+i+1<<"   "<<count+i+2<<endl; //支化点原子与主链上原子的键角
                cout<<count+comb[n].bonelen-2+bac+b<<"    "<<"1"<<"   "<<count+i+1<<"   "<<count+comb[n].bonelen+bac+1<<"  "<<count+comb[n].bonelen+bac+2<<endl;
                for(int j=1;j<comb[n].branchlen[b]-1;j++)
                {
                    cout<<count+comb[n].bonelen+bac+b-1<<"    "<<"1"<<"   "<<count+comb[n].bonelen+bac+1<<"   "<<count+comb[n].bonelen+bac+2<<"  "<<count+comb[n].bonelen+bac+3<<endl;
                    bac++;
                }
                bac=bac+2;
                cout<<count+comb[n].bonelen-2+bac+b-1<<"  "<<"1"<<"   "<<count+i<<"  "<<count+i+1<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]+1<<endl;
                cout<<count+comb[n].bonelen-2+bac+b<<"  "<<"1"<<"   "<<count+i+2<<"  "<<count+i+1<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]+1<<endl;
                b++;
            }
            else
            {
                cout<<count+i<<"  "<<"1"<<"   "<<count+i<<"   "<<count+i+1<<"   "<<count+i+2<<endl;// 主链上键角
            }
        }
        count+=comb[n].bonelen+bac;
    }

    count=0;
    mol=0;
    cout<<endl<<"Dihedrals"<<endl<<endl;  //二面角/////
    for(int n=0;n<chain;n++)
    {
        bac=0,b=1;
        for(int i=1;i<comb[n].bonelen-2;i++)
        {
            if(i+1==comb[n].branchpos[b])
            {
                cout<<count+i<<"  "<<"1"<<"   "<<count+i+mol<<"   "<<count+i+1+mol<<"   "<<count+i+2+mol<<"   "<<count+i+3+mol<<endl; //支化点原子与主链上原子的键角
                cout<<count+comb[n].bonelen-3+bac+b<<"    "<<"1"<<"   "<<count+i+1+mol<<"   "<<count+comb[n].bonelen+bac-b+2+mol<<"  "<<count+comb[n].bonelen+bac-b+3+mol<<"  "<<count+comb[n].bonelen+bac-b+4+mol<<endl;
                int bb=bac;
                for(int j=1;j<comb[n].branchlen[b]-2;j++)
                {
                    bac++;
                    cout<<count+comb[n].bonelen-3+bac+b<<"    "<<"1"<<"   "<<count+comb[n].bonelen+bac-b+1+mol<<"   "<<count+comb[n].bonelen+bac-b+2+mol<<"   "<<count+comb[n].bonelen+bac-b+3+mol<<"  "<<count+comb[n].bonelen+bac-b+4+mol<<endl;
                    
                }
                bac=bac+2;
                cout<<count+comb[n].bonelen-3+bac+b-1<<"  "<<"1"<<"   "<<count+i-1+mol<<"   "<<count+i+mol<<"  "<<count+i+1+mol<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]+3-b+mol<<endl;
                cout<<count+comb[n].bonelen-3+bac+b<<"  "<<"1"<<"   "<<count+i+3+mol<<"   "<<count+i+2+mol<<"  "<<count+i+1+mol<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]+3-b+mol<<endl;
                bac=bac+2;
                cout<<count+comb[n].bonelen-3+bac+b-1<<"  "<<"1"<<"   "<<count+i+mol<<"  "<<count+i+1+mol<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]-b+1+mol<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]-b+2+mol<<endl;
                cout<<count+comb[n].bonelen-3+bac+b<<"  "<<"1"<<"   "<<count+i+2+mol<<"  "<<count+i+1+mol<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]-b+1+mol<<"   "<<count+comb[n].bonelen+bac-comb[n].branchlen[b]-b+2+mol<<endl;
                b++;
            }
            else
            {
                cout<<count+i<<"    "<<"1"<<"     "<<count+i+mol<<"     "<<count+i+1+mol<<"       "<<count+i+2+mol<<"     "<<count+i+3+mol<<endl;// 主链上键角
            }
        }
        if(comb[n].branchnum==0)
        {
            count+=comb[n].bonelen-1;
        }
        else
        {
            count+=comb[n].bonelen+bac-2;
        }
        mol++;
    }



return 0;
}
