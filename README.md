# PE_Chain_for_MD
Generate a input United Atom file for PE MD in LAMMPS
可生成自定义链长、链数的聚乙烯联合原子模型初始构象
# branched_PE
  用于生成支化聚乙烯链，相邻的支化位点距离建议大于4个C,过多支化可能还存在Bug(当然，把每条链的描述信息中支链数设置为0，也可生成线性的链)
## chain_length.txt
  >第1行忽略
  >第2行表示链的总数N
  >第3行到第n+2行的第一个数表示每条链的长度，第二位数表示支链个数（Nb），往后2Nb个数表示第i个支化位点及其支化长度
## branchedpe.cpp
   源码
## branchedpe.exe
   编译结果
## fcc_unwrapped_coordinate.dat
  非周期性边界条件的坐标文件信息，用以作为初始构象输入LAMMPS
# linear_PE
  用于生成线性聚乙烯链，简化输入参数