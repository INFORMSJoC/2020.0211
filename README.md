[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Disjoint Bilinear Optimization: A Two-stage Robust Optimization Perspective


## Cite

To cite this material, please cite this repository, using the following DOI.

[![DOI](https://zenodo.org/badge/428188687.svg)](https://zenodo.org/badge/latestdoi/428188687)

Below is the BibTex for citing this version of the code.

```
@article{DBOTSROA21,
  author =        {N. Graham, H. Hu, J. Im, X. Li and H. Wolkowicz},
  publisher =     {INFORMS Journal on Computing},
  title =         {Disjoint Bilinear Optimization: A Two-stage Robust Optimization Perspective},
  year =          {2021},
  doi =           {10.5281/zenodo.5752594},
  url =           {https://github.com/INFORMSJoC/2020.0221}
}  
```

## Table of content
* [Requirements](#requirements)
* [Instances](#instances)
* [Algorithms](#Algorithms)

## Requirements
For this project, we use
* Matlab (the codes are written in Matlab R2020b)
* Yalmip (https://yalmip.github.io/) to pass optimization to solvers
* Gurobi (https://www.gurobi.com/) as the solver to solve Mixed Integer Linear Optimization problems
* Mosek (https://www.mosek.com/) to solve the conic optimization problems 

For licensing, see the LICENSE file. We emphasize that the user should respect the license of all the used sotfware, solvers, and packages. 

## Instances
This project contains four types of instances: Bimatrix, Constrained Bimatrix, Random Bilinear, and Convex Max. We explicitly mention the link between the data and the mathematical formulation for each class of instances here. Therefore, you can use Yalmip (https://yalmip.github.io/) and an appropriate solver to solve the instances.

To use R4B, you need to reformulate instances to a bilinear problem. In the paper, we explain how to do so for the constrained bimatrix games and the convex maximization problems. Then, you can directly call ``R4B`` function.

### Bimatrix
There are 50 randomly generated bimatrix games in this folder. Each of the instances contains *A_x* and *A_y*, which are the  payoff matrices of players, and *x_LH* and *y_LH*, which is a nash equiliburium. For the formulation of the optimization problem, see Problem (10) in the paper. 

### Constrained Bimatrix
This folder contains 50 randomly generated constrained bimatrix games. Each of the instances contains *A_x* and *A_y*, which are the  payoff matrices of players, *a_x*, *a_y*, *d_x*, and *d_y* to have the restrictions <img src="https://render.githubusercontent.com/render/math?math=a_x^T x \leq d_x"> and <img src="https://render.githubusercontent.com/render/math?math=a_y^T y \leq d_y">. The instance also contains *sol_x* and *sol_y*, which are the solutions for *x* and *y*, *sol_b_x* and *sol_b_y*, which are the solutions for *b_x* and *b_y*, and *nu_x* and *nu_y*, which are optimal values for  <img src="https://render.githubusercontent.com/render/math?math=\nu_x"> and <img src="https://render.githubusercontent.com/render/math?math=\nu_y">.

For the mathematical formulation of constrained bimatrix games, see Problem (11) in the paper.  

### Random BP
This folder contains 320 randomly generated blinear optimization instances. The name of each instance shows its characteristics. We consider 3 characteristics in a bilinear optimization problem: (i) number of the variables (which are denoted by <img src="https://render.githubusercontent.com/render/math?math=n_x"> for the vector variable *x* and by <img src="https://render.githubusercontent.com/render/math?math=n_y"> for the vector variable *y*); (ii) number of constraints (which are denoted by <img src="https://render.githubusercontent.com/render/math?math=m_x"> for the number of constraints on *x* and by <img src="https://render.githubusercontent.com/render/math?math=m_y"> for the number of constraints on *y*); (iii) the density of the matrices (denoted by <img src="https://render.githubusercontent.com/render/math?math=D_Q"> for the density of the matrix Q defining the objective fucntion, and <img src="https://render.githubusercontent.com/render/math?math=D_A"> for the density of the matrices defining the constraints). So, for example, the name BiL_mx20nx100my20ny100D_Q30D_A30Num1.mat indicates that this instance is for a bilinear optimization problem where dimensions of *x* and *y* are 100, and we have 20 constrains for *x* and 20 constraints for *y*, the density of the coefficient matrices for *x* and *y* is 0.3, and this is the first out of 5 instances with these charactestics.  

Each instance contains *Ax*, *Ay*, *bx*, *by*, and *Q*, to define the instance (see Problem (BO) in the paper). Also, the file contains *sol_x* and *sol_y*, which are the obtained optimal solution. We emphasize that the solutions are obtained with the relative optimality gap of 0.0001.

### Convex Max
This folder contains 13 instances for a convex maximization problem (see Problem (15) in the paper for the formultion). Each instance has the data for both the feasible region <img src="https://render.githubusercontent.com/render/math?math=\mathcal{X}_1"> and <img src="https://render.githubusercontent.com/render/math?math=\mathcal{X}_2">. In each file, everything with ubx1 comes from solving the convex relaxation with <img src="https://render.githubusercontent.com/render/math?math=\mathcal{X} = \mathcal{X}_1">. In this case, *xvalue_ubx1*, *yvalue_ubx1* and *Wvalue_ubx1* are therefore the values for *x*, *y* and *W* obtained from solving the convex relaxation with <img src="https://render.githubusercontent.com/render/math?math=\mathcal{X} = \mathcal{X}_1">. Similarly, *xvalue_ubx2*, *yvalue_ubx2* and *Wvalue_ubx2* are by considering <img src="https://render.githubusercontent.com/render/math?math=\mathcal{X} = \mathcal{X}_2">. 

*xvalue_lbx1*, *yvalue_lbx1* are then the solutions for *x* and *y* coming from solving the mountain climbing procedure in R4B. In the examples, *Jvert = [n_y]*. The data on *Q* and *c* are used to define the objective function, each row of the matrix *I* contains the indices in the set <img src="https://render.githubusercontent.com/render/math?math=\mathcal{J}_k">.


## Algorithms
We have implemeted the code for two classes of problems: bilinear optimization problems and convex maximization problems. 

For bilinear optimization problem, we have the function R4B inside the file R4B_linear. This function has the following form: 

``
[val,time,sol_x,sol_y,obj] = R4B(Q,A_x,b_x,A_y,b_y,tolerr,time_limit)
``

where *Q*, *A_x*, *b_x*, *A_y*, are *b_y* used in the definition of the optimization problem, *tolerr* is the relative optimality gap, and *time_limit* is the desired time limit. This function returns *val*, which is the obtained lowerbound, *time*, which is the solution time, *sol_x* and *sol_y* are the obtained solution, and *obj* is the obtained objective value. 

In the implemetation of this fucntion, we use *M_value* function to find the value of the M used in the reformulation (see PRoblem (MBC) in the paper). So, if the problem has a specific character, M can be obtained in a better way. In that case, we recommend the user to change this function. The proposed function is using the robust reformulation in Yalmip, which is slow.

The file, R4B_convexMAX.m also contains all the steps to construct the lowerbounds and upperbounds for the convex maximization problems. 

