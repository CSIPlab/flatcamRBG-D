# imageDepthLensless

This repository is the implementation of our paper: Joint Image and Depth Estimation with Mask-Based Lensless Cameras (IEEE transactions on computational imaging). The paper can be found on IEEE Xplore: https://ieeexplore.ieee.org/document/9144433.

## Abstract
![image](https://github.com/CSIPlab/imageDepthLensless/blob/master/doc/intro.png)
In this project, we propose an alternating gradient descent algorithm that jointly estimates a continuousdepth map and light distribution of the unknown scene from itslensless  measurements. We built a prototype lensless camera and present experimental results for reconstruction of intensity and depth maps of different real scenes. 

## Run the code
Run `demo_depth.m` file in Matlab.

## Parameters
Four example data files are included in this repository, they include a single card scene, two cards scene, a cup scene and a hand statue scene. The algorithm alternatively solves for image and depth map within a couple of loops. Here we list parameters that you can play with in the code.

- `dataset`: data files to be evaluated.
- `niter`: number of iterations in the loop.
- `damp`: l2 regularization parameter in the image solver.
- `itnlim`: number of iterations in the image solver.
- `maxIter`: number of iterations in the depth solver.
- `lambda`: TVL2 regularization parameter in the depth solver. 
