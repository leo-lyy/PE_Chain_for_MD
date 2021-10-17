（弱鸡的第一个repository）
# PE_Chain_for_MD
Generate a input molecular data file for MD in LAMMPS
可生成自定义链长、链数的聚乙烯初始构象
# 文件描述
## chain_length.txt
  >第1行忽略（来源于生成符合一定分子链长度分布的序列时的附加信息，在此处可忽略）
  >第2行表示链的数量n
  >第3行到第n+2表示每条链的长度
## PE_generator.cpp
   源码，第67行,cell这一变量用于设定盒子的大小，理论上建议设置为4的倍数。
## PE_generator.exe
   编译结果
## wrapped_coordinate.dat
  周期性边界条件的坐标文件信息
## unwrapped_coordinate.dat
  非周期性边界条件的坐标文件信息，用以作为初始构象输入LAMMPS
