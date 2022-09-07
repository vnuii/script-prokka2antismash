# Prokka to Antismash

## Introduction

该脚本前后共调用了prokka、convert2bed、mafft、hmmer、orthofisher、bedtools和antismash七个软件，可以将二代测序数据（细菌）通过prokka进行注释，用hmmer和目标蛋白进行比对，最后再用antismash进行目标次级代谢产物合成基因簇的分析。

## Flow Chart

![test](https://user-images.githubusercontent.com/103515667/188831579-e62e045c-424b-47f8-a37d-b5d3d467f37d.jpg)

## **!! Notice !!**

1. 运行该脚本需要两个数据，一个是细菌（最好是细菌吧）的**二代测序数据**，另一个是你需要比对的蛋白的**氨基酸序列**
1. 对于数据的格式有些许要求

​		**测序数据:** 换行符应为**Unix格式**而不是Windows格式，每个scallfold spades与>间无空格，						   如下

```shell
>scaffold_spades_01 #✅
> scaffold_spades_01 #❎
```

​		**蛋白数据:** 同上

​		两种文件的**后缀应都为.fasta**

3. 运行前需更改脚本中第15到17行的文件所在路径，具体👀⬇️

```shell
basepath_genomeseq=/home/von/ceshi/genomeseq #基因组测序数据所在的路径，注意使用绝对路径（下同）
base_path=/home/von/ceshi #你所希望的该脚本输出结果的储存路径，最好该路径包含基因组文件和目标蛋白文件（我懒得再看了，不包含应该也没事）
basepath_targetprtseq=/home/von/ceshi/targetprtseq #目标蛋白文件，即你后续准备blast的蛋白的氨基酸序列
```

4. 运行方式

```shell
nohup bash PATH/TO/prokka2antismash.sh &
exit #挺重要的步骤，如果你和我的情况一样
```

**这个是我本人的运行方式，我本人的情况是，登录服务器后，几分钟不动就会自动断开与服务器的连接。如果你和我的情况一样，可以采取上述步骤。**
