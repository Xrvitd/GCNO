# RFEPS: Reconstructing Feature-line Equipped Polygonal Surface 
Code of RFEPS.

Abs: Feature lines are important geometric cues in characterizing the structure of a CAD model. Despite great progress in both explicit reconstruction and implicit reconstruction, it remains a challenging task to reconstruct a polygonal surface equipped with feature lines, especially when the input point cloud is noisy and lacks faithful normal vectors. In this paper, we develop a multistage algorithm, named RFEPS, to address this challenge. The key steps include (1)denoising the point cloud based on the assumption of local planarity, (2)identifying the feature-line zone by optimization of discrete optimal transport, (3)augmenting the point set so that sufficiently many additional points are generated on potential geometry edges, and (4) generating a polygonal surface that interpolates the augmented point set based on restricted power diagram. We demonstrate through extensive experiments that RFEPS, benefiting from the edge-point augmentation and the feature-preserving explicit reconstruction, outperforms state-of-the-art methods in terms of the reconstruction quality, especially in terms of the ability to reconstruct missing feature lines.

Paper link: https://arxiv.org/abs/2212.03600
Doi: https://dl.acm.org/doi/10.1145/3550454.3555443

### Dependence

- CGAL 
- Eigen3
- Boost

### Makefile builds (Linux, other Unixes, and Mac. But we recommend using Windows.)

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

The example is in 'MAIN'. Include RFEPS in your project when testing and using it.

All the files is in 'RFEPS\data'. 

Please open ``OPENMP`` in Visual Studio to get the best performance.

The Restricted Power Diagram(RPD) in this project is a version that we implemented to facilitate debugging. If you want to get the fastest running speed, please use:
https://github.com/basselin7u/GPU-Restricted-Power-Diagrams


## Testing Platform
- Windows 10 
- Visual Studio 2022
- AMD Ryzen 5950X
- 64GB Momery

Thanks for BGAL library: https://github.com/BKHao/BGAL

If you use our code, please consider citing our work:
```
@article{Rui2022RFEPS,
  author    = {Xu, Rui and Wang, Zixiong and Dou, Zhiyang and Zong, Chen and Xin, Shiqing and Jiang, Mingyan and Ju, Tao and Tu, Changhe},
  title     = {RFEPS: Reconstructing Feature-line Equipped Polygonal Surface},
  journal   = {ACM Transactions on Graphics (TOG)},
  year      = {2022},
  publisher = {ACM},
  doi = {10.1145/3550454.3555443},
  booktitle = {ACM SIGGRAPH Asia 2022 Papers},
  numpages = {15},
  series = {SIGGRAPH Asia '22}
}
```

