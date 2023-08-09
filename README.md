# Network-based inference

## SDTNBI方法

子结构-药物-靶标网络推理方法（Substructure-Drug-Target Network based Inference, SDTNBI）是一种基于配体和网络的药物发现方法。SDTNBI可以在不依赖靶点的结构信息的情况下，通过来源于公共活性数据库、文献或专利中的药物-靶标相互作用对，构建相应网络并预测新药物的靶点。

首先，在开始阶段，资源传播次数为0（即：$`k = 0`$）。与$`D_i`$相连的所有邻居节点，也即$`D_i`$的各个己知靶标和己知子结构，会被各自分配到一个单位的初始资源。然后，在第一轮资源扩散过程中，每一个拥有初始资源的靶标节点，都会把它自己的资源平均地分给与之相连接的各个邻居节点。这时（即：$`k = 1`$时）可以发现，资源都移动到了药物节点上。在第二轮资源扩散过程中，每一个药物节点，又会把它自己的资源平均地分给与之相连的各个邻居节点。这时（即：$`k = 2`$时）又可以发现，资源又重新回到了靶标节点和子结构节点上。此时，对网络中的任意靶标（表示为$`T_j`$）而言，其拥有资源的多少，即是SDTNBI方法对药物一靶标相互作用$`D_i - T_j`$的预测打分。这一打分越高，表示$`D_j`$越有可能作用于$`T_j`$。通过这种方式，我们可以为网络中的所有药物进行系统性的药物一靶标相互作用预测。

Zengrui Wu and others, SDTNBI: an integrated network and chemoinformatics tool for systematic prediction of drug–target interactions and drug repositioning, Briefings in Bioinformatics, Volume 18, Issue 2, March 2017, Pages 333–347, https://doi.org/10.1093/bib/bbw012

## 安装依赖

* python >= 3.6
* numpy >= 1.16.2
* scipy >= 1.2.1
* rdkit >= 20.19.1.23

## SDTNBI使用教程

### 使用接口预测


### 载入默认网络并预测

```python
from network.utils import load_network, load_pretrained_network

smis = ["CN(C)C(=O)c1ccccc1C(=O)N(C)C",
    "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1",
    "CCCCNCc1ccc(S(=O)(=O)c2ccc(S(N)(=O)=O)s2)cc1",
    "COc1c(Cl)cc(Cl)c(OC)c1O",
    "Cc1cc(=O)[nH]c(=O)[nH]1",
    "O=c1[nH]c(=O)n(Cc2ccc(Cl)cc2)cc1F",
    "COc1ncnc2nccnc12",
    "CCC(CC)(Cc1ccc(C)cc1)C(=O)NO",
    "COc1ccc2cc(C(C)C(=O)OC3COC4C(O)COC34)ccc2c1",
    "CC(=O)Nc1cc(N[S+2]([O-])([O-])C(F)(F)F)c(C)cc1C",
    "O=C(O)c1cccnc1"]

model = load_pretrained_network("chembl") # load pretrained model directly
# model = load_network("data/chembl-sdtnbi.pkl.gz") # load your trained model
F = model.predict(chemicals=smis, target='Q99835')
print(F)
```

### 训练自己的网络

#### 以默认参数训练网络

```python
from network.sdtnbi import SDTNBI
from network.utils import save_network

model = SDTNBI()
model.train(network="./data/drugbank-network.csv", drugs="./data/drugbank-drug.csv")

smis = ["CN(C)C(=O)c1ccccc1C(=O)N(C)C",
    "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1",
    "CCCCNCc1ccc(S(=O)(=O)c2ccc(S(N)(=O)=O)s2)cc1",
    "COc1c(Cl)cc(Cl)c(OC)c1O",
    "Cc1cc(=O)[nH]c(=O)[nH]1",
    "O=c1[nH]c(=O)n(Cc2ccc(Cl)cc2)cc1F",
    "COc1ncnc2nccnc12",
    "CCC(CC)(Cc1ccc(C)cc1)C(=O)NO",
    "COc1ccc2cc(C(C)C(=O)OC3COC4C(O)COC34)ccc2c1",
    "CC(=O)Nc1cc(N[S+2]([O-])([O-])C(F)(F)F)c(C)cc1C",
    "O=C(O)c1cccnc1"]

F = model.predict(chemicals=smis, target='Q99835')
# save_network(model, "data/small-sdtnbi.pkl.gz") # uncomma for saving model

print(F)
```

#### 更改网络子结构参数并训练网络

当训练较小的网络时，不需要使用分子指纹（子结构数）默认参数（3，16384），可以使用较小的分子指纹（2，2048）

```python
from sdtnbi import *

from network.sdtnbi import SDTNBI
from network.utils import save_network

model = SDTNBI()
model.train(network="./data/drugbank-network.csv",
    drugs="./data/drugbank-drug.csv", r=2, nbits=2048)

smis = ["CN(C)C(=O)c1ccccc1C(=O)N(C)C",
    "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1",
    "CCCCNCc1ccc(S(=O)(=O)c2ccc(S(N)(=O)=O)s2)cc1",
    "COc1c(Cl)cc(Cl)c(OC)c1O",
    "Cc1cc(=O)[nH]c(=O)[nH]1",
    "O=c1[nH]c(=O)n(Cc2ccc(Cl)cc2)cc1F",
    "COc1ncnc2nccnc12",
    "CCC(CC)(Cc1ccc(C)cc1)C(=O)NO",
    "COc1ccc2cc(C(C)C(=O)OC3COC4C(O)COC34)ccc2c1",
    "CC(=O)Nc1cc(N[S+2]([O-])([O-])C(F)(F)F)c(C)cc1C",
    "O=C(O)c1cccnc1"]

F = model.predict(chemicals=smis, target='Q99835')
# save_network(model, "data/small-sdtnbi.pkl.gz") # uncomma for saving model

print(F)
```

## NBI方法的矩阵描述

在介绍SDTNBI方法前，我们需要先了解NBI方法。我们首先需要定义2个集合：药物集合$`D = \{D_1, D_2, ..., D_{N_D}\}`$，靶标集合$`T = \{T_1, T_2, ..., T_{N_T}\}`$以及一个子结构集合$`S = \{S_1, S_2, ..., S_{N_S}\}`$。三者的关系可以表示为两个矩阵形式：

```math
A_{DT} = \begin{cases}
    1, & \text{if $D_i$ acts on $T_j$} \\
    0, & \text{otherwise}
\end{cases}
```

```math
A_{DS} = \begin{cases}
    1, & \text{if $D_i$ contains $S_j$} \\
    0, & \text{otherwise}
\end{cases}
```

基于以上矩阵，初始资源矩阵可以表示为：

```math
A = \begin{bmatrix}
    0 & A_{DS} & A_{DT} \\
    A_{DS}^T & 0 & 0 \\
    A_{DT}^T & 0 & 0 \\
\end{bmatrix}
```

扩散矩阵W可以通过以下公式计算：

```math
W_{i, j} = \begin{cases}
    \frac{A(i, j)}{\sum_{l=1}^{N_S + N_D + N_T}A(i, l)} & \text{if $A(i, j) \neq 0$} \\
    0, & \text{otherwise}
\end{cases}
```

假设资源扩散的次数是$k$，那么最终的资源矩阵就可以通过如下公式计算：

```math
F = A \times W^k
```

最终矩阵中，$`F(i, N_D + N_S + j)`$的值即使SDTNBI方法给出的药物-靶标相互作用$`D_i - T_j`$的预测打分，其中$`i \in (0, N_D]`$、$`j \in (0, N_T]`$

## SDTNBI方法的矩阵描述

对于SDTNBI方法，我们需要在原有集合的基础上定义4个集合：药物集合$`D = \{D_1, D_2, ..., D_{N_D}\}`$，靶标集合$`T = \{T_1, T_2, ..., T_{N_T}\}`$，一个子结构集合$`S = \{S_1, S_2, ..., S_{N_S}\}`$以及未知活性化合物集合$`C = \{C_1, C_2, ..., C_{N_C}\}`$。三者的关系可以表示为三个矩阵形式：

```math
A_{DT} = \begin{cases}
    1, & \text{if $D_i$ acts on $T_j$} \\
    0, & \text{otherwise}
\end{cases}
```

```math
A_{DS} = \begin{cases}
    1, & \text{if $D_i$ contains $S_j$} \\
    0, & \text{otherwise}
\end{cases}
```

```math
A_{CS} = \begin{cases}
    1, & \text{if $C_i$ contains $S_j$} \\
    0, & \text{otherwise}
\end{cases}
```

基于以上矩阵，初始资源矩阵可以表示为：

```math
A = \begin{bmatrix}
    0 & 0 & A_{CS} & 0 \\
    0 & 0 & A_{DS} & A_{DT} \\
    A_{CS}^T & A_{DS}^T & 0 & 0 \\
    0 & A_{DT}^T & 0 & 0 \\
\end{bmatrix}
```

我们可以进一步排除其中那些<font color=green>没有己知靶标信息的分子</font>的影响。获得仅有药物-靶标的资源矩阵：

```math
A\_ = \begin{bmatrix}
    0 & 0 & 0 & 0 \\
    0 & 0 & A_{DS} & A_{DT} \\
    0 & A_{DS}^T & 0 & 0 \\
    0 & A_{DT}^T & 0 & 0 \\
\end{bmatrix}
```

基于仅有药物-靶标的资源矩阵，我们仍可以计算扩散矩阵：

```math
W_{i, j} = \begin{cases}
    \frac{A\_(i, j)}{\sum_{l=1}^{N_S + N_D + N_T}A\_(i, l)} & \text{if $A\_(i, j) \neq 0$  $\sum_{l=1}^{N_S + N_D + N_T}A\_(i, l) \neq 0$} \\
    0,  \text{otherwise}
\end{cases}
```

假设资源扩散的次数是$`k`$，那么最终的资源矩阵就可以通过如下公式计算：

```math
F = A \times W^k
```

根据论文结果，我们取扩散次数为2

最终矩阵中，$`F(i, N_D + N_S + j)`$的值即使SDTNBI方法给出的药物-靶标相互作用$`D_i - T_j`$的预测打分，其中$`i \in (0, N_D]`$、$`j \in (0, N_T]`$

## SDTNBI矩阵运算的简化

我们尝试对矩阵进行拆解：

```math
\begin{aligned}
  F &=  A \times W^2\\
    &=  \begin{bmatrix}
        0 & 0 & A_{CS} & 0 \\
        0 & 0 & A_{DS} & A_{DT} \\
        A_{CS}^T & A_{DS}^T & 0 & 0 \\
        0 & A_{DT}^T & 0 & 0 \\
        \end{bmatrix}
        \times
        \begin{bmatrix}
        0 & 0 & 0 & 0 \\
        0 & 0 & W_{DS}^{\prime} & W_{DT}^{\prime} \\
        0 & W_{DS}^{\prime \prime} & 0 & 0 \\
        0 & W_{DT}^{\prime \prime} & 0 & 0 \\
        \end{bmatrix}
        \times
        \begin{bmatrix}
        0 & 0 & 0 & 0 \\
        0 & 0 & W_{DS}^{\prime} & W_{DT}^{\prime} \\
        0 & W_{DS}^{\prime \prime} & 0 & 0 \\
        0 & W_{DT}^{\prime \prime} & 0 & 0 \\
        \end{bmatrix}\\  
\end{aligned}
```

我们只关心右上角DT矩阵的计算结果：

```math
\begin{aligned}
  F_{DT} &= \begin{bmatrix}
            0 & 0 & A_{CS} & 0 \\
            \end{bmatrix}
            \times
            \begin{bmatrix}
            0 & 0 & 0 & 0 \\
            0 & 0 & W_{DS}^{\prime} & W_{DT}^{\prime} \\
            0 & W_{DS}^{\prime \prime} & 0 & 0 \\
            0 & W_{DT}^{\prime \prime} & 0 & 0 \\
            \end{bmatrix}
            \times
            \begin{bmatrix}
            0 & W_{DT}^{\prime T} & 0 & 0\\
            \end{bmatrix}^T\\
         &= \begin{bmatrix}
            0 & A_{CS} \times W_{DS}^{\prime \prime T} & 0 & 0 \\
            \end{bmatrix}
            \times
            \begin{bmatrix}
            0 & W_{DT}^{\prime T} & 0 & 0 \\
            \end{bmatrix}^T\\
         &= A_{CS} \times W_{DS}^{\prime \prime T} \times W_{DT}^{\prime}
\end{aligned}
```

python实现如下：

```python
%%time
import numpy as np

N_C = 100
N_D = 400
N_S = 2048
N_T = 300
k = 2

A_CS = np.matrix(np.random.randn(N_C, N_S)+1., dtype=np.float16) # dim: 1x4

A_DT = np.matrix(np.random.randn(N_D, N_T)+1., dtype=np.float16) # dim: 4 x 5

A_DS = np.matrix(np.random.randn(N_D, N_S)+1., dtype=np.float16) # dim: 4 x 4



W_DS_T = A_DS.T/np.sum(A_DS, axis=0).T

W_DT = A_DT/(np.sum(A_DT, axis=1)+np.sum(A_DS, axis=1))

F = np.dot(np.dot(A_CS, W_DS_T), W_DT)
```

## ChemBL数据库清洗流程

* 首先使用SQL进行初步筛选，只保留活性（IC50,EC50,Ki,Kd）小于1umol的化合物进行进一步筛选
* 仅保留靶点为人源单一蛋白的活性数据
* 性质筛选：-0.4<=alogP<=5.6, rtb<=10, HBA<=10, HBD<=5, 180<=MW<=500
* 来源为文献数据，有doi号

最终筛选得到了154529个有效活性数据对，其中活性化合物个140641，靶标1211个。

进一步使用RDkit进行聚类分析，对每个靶点的化合物进行以相似性cut-off为0.6，基于MorganFP的Butina聚类，最后保留各个聚类的中心分子，最终减少到36010个有效数据对，其中活性化合物个35428个，靶标1211个。
