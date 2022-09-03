# RFEPS: Reconstructing Feature-line Equipped Polygonal Surface (Working, temporarily unavailable)
Code of RFEPS (SIGGRAPH Asia 2022 Conditionally Accepted).

Thanks for the simple and easy to use BGAL library: https://github.com/BKHao/BGAL

### Dependence

- CGAL 
- Eigen3
- Boost

### Makefile builds (Linux, other Unixes, and Mac)

```
git clone https://github.com/Xrvitd/RFEPS
cd RFEPS
mkdir build && cd build
cmake ..
make -j8
make install
```


### MSVC on Windows

```
git clone https://github.com/Xrvitd/RFEPS
```
Open cmake-gui

```
Where is the source code: RFEPS

Where to build the binaries: RFEPS/build
```

note: check the location of dependencies and install. It is recommended to use vcpkg to add dependencies.

Configure->Generate->Open Project

ALL_BUILD->INSTALL

## Test

The example and data is in 'test'. Include RFEPS in your project when testing and using it.

The Restricted Power Diagram(RPD) in this project is a version that we implemented to facilitate debugging. If you want to get the fastest running speed, please use:
https://github.com/basselin7u/GPU-Restricted-Power-Diagrams



If you use our code, please consider citing our work:
```
@article{Rui2022RFEPS,
  author    = {Rui Xu, Zixiong Wang, Zhiyang Dou, Zong Chen, Shiqing Xin, Mingyan Jiang, Tao Ju, Changhe Tu},
  title     = {RFEPS: Reconstructing Feature-line Equipped Polygonal Surface},
  journal   = {ACM Transactions on Graphics (TOG)},
  year      = {2022},
  publisher = {ACM}
}
```
Coming soon.
