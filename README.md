# Globally Consistent Normal Orientation for Point Clouds by Regularizing the Winding-Number Field

This paper is published at ACM Transactions on Graphics (SIGGRAPH 2023). Please refer to our [paper](https://arxiv.org/abs/2304.11605) and [project](https://xrvitd.github.io/Projects/GCNO/index.html) for more detail.

We just got the SIGGRAPH 2023 Best Paper Award!!! Congrats!

![](./pics/teaser4.png)


```
@article{Rui2023GCNO,
author = {Xu, Rui and Dou, Zhiyang and Wang, Ningna and Xin, Shiqing and Chen, Shuangmin and Jiang, Mingyan and Guo, Xiaohu and Wang, Wenping and Tu, Changhe},
title = {Globally Consistent Normal Orientation for Point Clouds by Regularizing the Winding-Number Field},
year = {2023},
issue_date = {August 2023},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {42},
number = {4},
issn = {0730-0301},
url = {https://doi.org/10.1145/3592129},
doi = {10.1145/3592129},
month = {jul},
articleno = {111},
numpages = {15},
keywords = {Voronoi diagram, winding number, normal orientation, raw point cloud, optimization}
}
```

## Tested Platform

- Windows 10 
- Visual Studio 2022 Professional
- AMD Ryzen 5950X
- 32GB Memory

## Dependencies
- CGAL 
- Eigen3
- Boost


**Important: Please using vcpkg to install dependent libraries! And please use  "git clone" to install vcpkg, otherwise you may get errors during CGAL installation.**

- .\vcpkg install boost:x64-windows

- .\vcpkg install cgal:x64-windows

  ​	use "git pull" if you get errors with the "gmp" library.

- .\vcpkg install Eigen3:x64-windows

- .\vcpkg integrate install



## MSVC on Windows

```
Download this project: NormalOrientation
```
Open Cmake-GUI

```
Where is the source code: NormalOrientation

Where to build the binaries: NormalOrientation/build
```

**Note: check the location of dependencies and install. It is recommended to use vcpkg to add dependencies.**

```
Configure->Generate->Open Project

ALL_BUILD -> Build
Turn Debug to Release -> ALL_BUILD -> Build
```

Please set `MAIN` as Startup Project, and make the following changes:

```
Properties -> Configuration Properties -> C/C++ -> Code Generation -> 
- Enable Parallel Code Generation : Yes
- Enable Enhanced Instruction Set : AVX2
- Floating Point Model : Fast

Properties -> Configuration Properties -> C/C++ -> Language -> Open MP Support : Yes
```

## Test

- All examples are in `MAIN`. 
- All the files is in `NormalOrientation\data`. 
- The output files is in `NormalOrientation\data\Out\`.
- The result of operations are in `NormalOrientation\data\MyResult\`, for checking whether the program is running correctly.


### Important Tips: 

💡💡💡 **Speed**
  
This code is not optimized for speed, but for clarity. Please open Openmp and AVX2 in Visual Studio to speed up the code.
Please set the floating point model to fast in Visual Studio to speed up the code.
The default number of Openmp parallel threads is 28, set according to an AMD Ryzen 5950x CPU, 
please set different number of threads according to the CPU you use to get the best running effect.

<img src="pics\image-20230116154151378.png" alt="image-20230116154151378" style="zoom:40%;" />

<img src="pics\image-20230116154210583.png" alt="image-20230116154210583" style="zoom:40%;" />



💡💡💡 **Stop condition**

For viewing the optimization process in more detail, we have not set the optimization stop condition. Please manually stop the optimization and view all iteration results in the `data\out` folder


## License
GCNO is under [AGPL-3.0](https://www.gnu.org/licenses/agpl-3.0.en.html), so any downstream solution and products (including cloud services) that include GCNO code inside it should be open-sourced to comply with the AGPL conditions. For learning purposes only and not for commercial use. If you want to use it for commercial purposes, please contact us first.

