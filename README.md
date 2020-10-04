# Joint Image and Depth Estimation with Mask-Based Lensless Cameras

This repository is the implementation of our paper: [Joint Image and Depth Estimation with Mask-Based Lensless Cameras](https://ieeexplore.ieee.org/document/9144433).

Please cite the following paper when using this code or data:
```
@ARTICLE{9144433,
  author={Y. {Zheng} and M. {Salman Asif}},
  journal={IEEE Transactions on Computational Imaging}, 
  title={Joint Image and Depth Estimation With Mask-Based Lensless Cameras}, 
  year={2020},
  volume={6},
  number={},
  pages={1167-1178},}
```

## Abstract
![image](https://github.com/CSIPlab/imageDepthLensless/blob/master/doc/intro.png)

In this project, we propose an alternating gradient descent algorithm that jointly estimates a continuousdepth map and light distribution of the unknown scene from itslensless  measurements. We built a prototype lensless camera and present experimental results for reconstruction of intensity and depth maps of different real scenes. 

## Run the code
Run `demo_depth.m` file in Matlab.

## Expected results
The algorithm has 3 main steps, here we list the expected results in each step.
### Step 1: Initialization
![image](https://github.com/CSIPlab/imageDepthLensless/blob/master/doc/depth_sweep.gif)

In this step, we try to reconstruct the scene at 10 predefined depth planes separately, then we pick the one that has the closest approximation as our initialization.

### Step 2: Optimization
![image](https://github.com/CSIPlab/imageDepthLensless/blob/master/doc/depth_est.gif)

In this step, we jointly estimate image and depth. In each loop, we will first estimate the depth value of each pixel using gradient descent based on the current estimated image, then estimate the image intensity based on the estimated depth.

### Step 3: Final results
![image](https://github.com/CSIPlab/imageDepthLensless/blob/master/doc/imgDepth.png)

In this step, we cast away the depth values that belong to pixels with low image intensities, because they it is inacurate to estimate the depth of the low light components using our method.

## Parameters
Four example data files are included in this repository, they include a single card scene, two cards scene, a cup scene and a hand statue scene. The algorithm alternatively solves for image and depth map within a couple of loops. Here we list some parameters that you can play with in the code.

- `dataset`: data files to be evaluated.
- `niter`: number of iterations in the loop.
- `damp`: l2 regularization parameter in the image solver.
- `itnlim`: number of iterations in the image solver.
- `maxIter`: number of iterations in the depth solver.
- `lambda`: TVL2 regularization parameter in the depth solver. 

## Packages
We include the packages we used in the /code/utils/ folder. Here we list these packages and their sources.

- minFunc: https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
- LSQR: https://web.stanford.edu/group/SOL/software/lsqr/
- export_fig: https://github.com/altmany/export_fig
