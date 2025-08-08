<div align=left><img src="PHengLEI\Documents/PHengLEI-LOGO.png" style="zoom: 15%;"/><br>

# 风雷软件（PHengLEI2506.v1818）

备注：主分支为master分支，更稳定，更多新功能请切换到ActiveBranch分支。

## 1.软件简介

&ensp;&ensp;&ensp;&ensp;风雷软件<sup>[1,2]</sup>（PHengLEI，Platform for Hybrid ENGineering simulation of flows）是中国空气动力研究与发展中心（CARDC）研发的面向流体工程的混合CFD平台。平台以面向对象的设计理念，采用C++语言编程。2020年12月，风雷软件正式面向全国开源，与其他开源CFD软件相比，风雷软件具有扩展能力强、开发难度低、计算效率高等特点。（更多介绍请阅读PHengLEI/Documents文件夹下的《风雷软件应用与开发指南》）

&ensp;&ensp;&ensp;&ensp;风雷软件更多动态和Demo请登录以下网址查看：

- 官网地址: https://www.cardc.cn/nnw/products.aspx?t=9
- 代码库地址：https://forge.osredm.com/PHengLEI/PHengLEI
- 算例库地址：https://forge.osredm.com/PHengLEI/PHengLEI-TestCases
- 论坛（常见问题、算例展示、技术分享）地址：https://osredm.com/forums/theme/38
- 风雷软件用户录制视频教程地址：https://www.bilibili.com/video/BV1eX4y1T7yW
[1] 赵钟,等.风雷（PHengLEI）通用CFD软件设计[J]. 计算机工程与科学, 2020, 42(2): 210-219.( ZHAO Zhong, et al. Design of general CFD software PHengLEI [J]. Computer Engineering & Science, 2020, 42(2): 210-219. (in Chinese) )
[2] 赵钟,等.适用于任意网格的大规模并行CFD计算框架PHengLEI[J]. 计算机学报, 2019, 42(11):2368
-2383. ( ZHAO Zhong, et al. PHengLEI: A Large Scale Parallel CFD Framework for Arbitrary Grids [J]. Chinese Journal of Computers, 2019, 42(11): 2368-2383. (in Chinese) )

**声明：** 若用户将该软件用于学术研究或工程应用，须在相关的论文成果的显要位置处标注基于“风雷（PHengLEI）“软件，并引用“风雷（PHengLEI）”软件相关的参考文献（例如[1]和[2]）。

## 2.软件功能

&ensp;&ensp;&ensp;&ensp;风雷软件是一款结构/非结构通用CFD软件，计算范围覆盖低速、亚跨和高速；发展面向多区域多算法融合的混合计算框架，支持结构/非结构不同网格类型混合，二阶/高阶不同求解精度混合，可压/不可压不同求解器类型混合等多种混合策略；构建适配任意硬件架构的多层级混合并行模型，支持MPI、MPI/OpenMP、MPI/CUDA等多种同构及异构并行策略，具备百亿网格百万核并行模拟能力；发展多种高效计算加速策略，支持LUSGS/GMRES等多种高效隐式加速算法、具备湍流AI预测能力；支持高精度方法、SIMPLE不可压、多组分计算、多孔介质模拟，提供SA/SST/湍流壁面函数等多种湍流模型、具备RANS/DES/LES的定常/非定常湍流模拟能力，可满足常规工程应用需求；支持网格加密、网格变形、重叠网格、动网格等多种网格功能；支持LBM、基于LBM的浸入边界流固耦合方法，可满足用户多样化需求。同时，风雷软件也提供常用前、后置软件接口，如Gridgen、Pointwise、ICEM-CFD、FieldView、Tecplot、Ensight、ParaView等。


## 3.代码获取

1. 环境准备，安装git，官网地址：https://git-scm.com/;
2. 进入风雷代码库，点击右上角Fork按钮；
3. Fork完成后，会生成并跳转到新的仓库，复制新仓库地址，如https://www.osredm.com/p68217053/PHengLEI;
4. 在本地选择一个目录，右键打开git bash，输入命令进行代码克隆，如git clone https://www.osredm.com/p68217053/PHengLEI;
5. 输入用户名和密码，其中用户名是指上面命令中p开头的用户名，比如p68217053；
6. 项目克隆完成后，进入项目目录，默认分支为master分支，可通过命令切换到ActiveBranch分支，如git checkout ActiveBranch;
7. 切换到开发分支后，可输入git log查看日志，确认当前版本。

## 4.软件安装

&ensp;&ensp;&ensp;&ensp;风雷软件能够在Windows、Linux、Mac系统下运行，源代码采用C++语言编写，需要CMake软件构建项目，并行计算采用MPI库。因此，操作系统必须提供C++编译器、CMake2.8以上版本软件和MPI1.0或MPI2.0标准库。

### 4.1 Windows环境配置

&ensp;&ensp;&ensp;&ensp;Windows环境下所有必备软件按照默认步骤安装即可。

1. 安装Microsoft Visual Studio 2012以上版本；
2. 安装MPI库，推荐采用MSMPI；
3. 安装Cmake。

### 4.2 Linux环境配置

#### 4.2.1 本地构建模式

1. 安装基础工具链：
•	GCC 编译器
•	MPICH3 库（MPI 实现）
•	CMake（最低版本要求 3.18）
2. 将上述工具的可执行路径添加至环境变量（如 ~/.bashrc 或~/.zshrc）
3. 在项目根目录 PHengLEI/ 下执行make config
此命令触发：
•	第三方依赖库的自动化编译
•	风雷求解器的一键构建
**注意：**
&ensp;&ensp;&ensp;&ensp;构建系统集成编译环境自检功能。若检测到非推荐环境（如低版本编译器），需手动确认（输入 y + 回车）方可继续。

#### 4.2.2 预配置超算构建模式

1. 在项目根目录 PHengLEI/ 下执行：
<BASH>
sh ./hpc_build.sh
2. 根据提示选择目标超算平台及编译环境，当前支持：

| 超算平台 | 编译环境 |
|:------------:|:------------:|
| 成都超算中心   | GCC/Intel   |
| 天津超算中心   | ARM       |
| 广州超算中心   | GCC   |
| 无锡超算中心   | GCC   |
| 济南超算中心   | Intel   |

**注意：**
&ensp;&ensp;&ensp;&ensp;构建完成后，必须执行以下命令激活环境配置：
&ensp;&ensp;&ensp;&ensp;<BASH>
&ensp;&ensp;&ensp;&ensp;source env.sh
&ensp;&ensp;&ensp;&ensp;否则后续作业提交可能因环境变量未加载而失败。

**附加说明**
•	本系统同时支持逐步手动编译流程，详见PHengLEI/Documents文件夹下的《风雷软件应用与开发指南》。
•	自动化构建基于 Makefile 统一构建系统，实现依赖解析、环境检测与编译流程集成。

## 5.如何贡献

&ensp;&ensp;&ensp;&ensp;针对基于风雷软件进行二次开发的用户，后续添加如何测试、提交Pull Request的步骤。

&ensp;&ensp;&ensp;&ensp;软件贡献者信息参见项目根目录下的Contributing文件（待增加）。

## 6.软件版权

&ensp;&ensp;&ensp;&ensp;风雷软件开源协议参见项目根目录下的LICENSE文件。

## 7.鸣谢

&ensp;&ensp;&ensp;&ensp;感谢所有参与风雷软件开发与推广的工作人员，也感谢所有支持风雷软件发展并提出宝贵意见和建议的广大用户。

## 8.联系我们

&ensp;&ensp;&ensp;&ensp;电子邮箱：phenglei@126.com。