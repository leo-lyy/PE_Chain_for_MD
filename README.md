（弱鸡的第一个repository）
# PE_Chain_for_MD
Generate a input molecular data file for MD in LAMMPS
可生成自定义链长、链数以及少量支化的聚乙烯初始构象

# Linear_pe 线性聚乙烯
此文件夹下代码用于生成线性分子链
## chain_length.txt
  >分子链特征描述文件
  >第1行忽略（来源于生成符合一定分子链长度分布的序列时的附加信息，在此处可忽略）
  >第2行表示链的数量n
  >第3行到第n+2表示每条链的长度
  
## fcc_pe.cpp
  > 源码
  
## fcc_pe.exe
  > 编译结果
  
## wrapped_coordinate.dat
  > 周期性边界条件的坐标文件信息
  
## unwrapped_coordinate.dat
  > 非周期性边界条件的坐标文件信息，用以作为初始构象输入LAMMPS

# branched_pe 支化聚乙烯
此文件夹下代码用于生成支化分子链
当前代码部分完成，能够生成每条分子链支链数量小于等于2
支链数量大于等于3时，二面角关系会出错
## chain_length.txt
  > 分子链特征描述文件
  > 第二行表示分子的数量 N
  > 下面 N 行中， 第一个数 a 是支链数量，后面开始第2*i-1个数表示支化位点，第2*i个数表示该位点上的支化长度

## branchedpe.cpp
  > 源码
  
## branchedpe.exe
  > 编译结果
  
