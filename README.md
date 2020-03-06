# RDDP_PowerSystem

This toolbox implemented the Robust Dual Dynamic Programming(RDDP) for solving real-time hydro-thermal dispatch problems under the wind power uncertainty.

The RDDP solves multi-stage minimax problems in the following form:

![equation](http://latex.codecogs.com/gif.latex?$\operatorname{minimize} \quad \boldsymbol{q}_{1}^{\top} \boldsymbol{x}_{1}+\mathcal{Q}_{2}\left(\boldsymbol{x}_{1}\right)$)
subject to ![equation](http://latex.codecogs.com/gif.latex?$\quad W_{1} x_{1} \geq h_{1}$)
![equation](http://latex.codecogs.com/gif.latex?$x_{1} \in \mathbb{R}^{n_{1}}$)

In which ![equation](http://latex.codecogs.com/gif.latex?$\mathcal{Q}_{t}\left(\boldsymbol{x}_{t-1}\right)$) is called the stage *t* worst-case cost-to-go function, defined as below:

![equation](http://latex.codecogs.com/gif.latex?$\mathcal{Q}_{t}\left(\boldsymbol{x}_{t-1}\right) = \max \left\{Q_{t}\left(\boldsymbol{x}_{t-1} = \boldsymbol{\xi}_{t}\right): \boldsymbol{\xi}_{t} \in \boldsymbol{\Xi}_{t}\right\}$)

![equation](http://latex.codecogs.com/gif.latex?$\begin{aligned}
Q_{t}\left(\boldsymbol{x}_{t-1} ; \boldsymbol{\xi}_{t}\right)=\left[\begin{array}{cc}
\operatorname{minimize} & \boldsymbol{q}_{t}^{\top} \boldsymbol{x}_{t}+\mathscr{Q}_{t+1}\left(\boldsymbol{x}_{t}\right) \\
\text { subject to } & \boldsymbol{T}_{t}\left(\boldsymbol{\xi}_{t}\right) \boldsymbol{x}_{t-1}+\boldsymbol{W}_{t} \boldsymbol{x}_{t} \geq \boldsymbol{h}_{t}\left(\boldsymbol{\xi}_{t}\right) \\
\boldsymbol{x}_{t} \in \mathbb{R}^{n_{t}}
\end{array}\right] \\
& t=2, \ldots, T
\end{aligned}$)

The worst-case cost-to-go function can be bounded from above and from below.
**Upper Bound**

![equation](http://latex.codecogs.com/gif.latex?$\bar{Q}_{t}\left(x_{t-1}\right)$
$=\operatorname{maximize}\left[\begin{array}{cc}\operatorname{minimize} & \boldsymbol{q}_{t}^{\top} \boldsymbol{x}_{t}+\bar{\alpha}_{t+1}\left(\boldsymbol{x}_{t}\right) \\ \text { subject to } & \boldsymbol{T}_{t}\left(\boldsymbol{\xi}_{t}\right) \boldsymbol{x}_{t-1}+\boldsymbol{W}_{t} \boldsymbol{x}_{t} \geq \boldsymbol{h}_{t}\left(\boldsymbol{\xi}_{t}\right) \\ & \boldsymbol{x}_{t} \in \mathbb{R}^{n_{t}}\end{array}\right]$
subject to $\xi_{t} \in \Xi_{t}$)

**Lower Bound**

![equation](http://latex.codecogs.com/gif.latex?$\underline{Q}_{t}\left(x_{t-1} ; \xi_{t}\right)=\left[\begin{array}{cl}
\operatorname{minimize} & \boldsymbol{q}_{t}^{\top} \boldsymbol{x}_{t}+\underline{Q}_{t+1}\left(\boldsymbol{x}_{t}\right) \\
\text { subject to } & \boldsymbol{T}_{t}\left(\boldsymbol{\xi}_{t}\right) \boldsymbol{x}_{t-1}+\boldsymbol{W}_{t} \boldsymbol{x}_{t} \geq \boldsymbol{h}_{t}\left(\boldsymbol{\xi}_{t}\right) \\
\boldsymbol{x}_{t} \in \mathbb{R}^{n_{t}}
\end{array}\right]$)



A test for modifyed IEEE 5-bus system(two wind farms and two hydro-units are added).
--------------------------------------------------------------------------------

                      RDDP_PS.jl (c) Shi Yunhui, 2020

 Iteration    UpperBound       LowerBound        Gap       Time (s)       # Solves
        1    5.675215e+03   5.675215e+03   0.000000e+00   2.570088e-01        912
        2    1.686230e+06   2.805209e+04   9.833640e-01   4.329488e-01       2698
        3    1.686230e+06   3.571683e+04   9.788185e-01   6.149280e-01       4484
        4    1.666846e+06   4.945067e+04   9.703328e-01   8.059549e-01       6270
        5    1.412543e+06   4.493406e+04   9.681892e-01   9.957180e-01       8056
        6    1.393159e+06   4.945067e+04   9.645047e-01   1.189456e+00       9842
        7    1.382754e+06   4.493406e+04   9.675039e-01   1.386206e+00      11628
        8    1.366194e+06   2.161969e+05   8.417524e-01   1.577504e+00      13414
        9    1.363603e+06   3.036965e+05   7.772837e-01   1.780142e+00      15200
       10    1.350337e+06   3.153418e+05   7.664718e-01   1.971623e+00      16986
       11    1.010938e+06   3.598488e+05   6.440446e-01   2.174286e+00      18772
       12    8.427517e+05   3.894124e+05   5.379275e-01   2.377907e+00      20558
       13    8.233682e+05   4.360300e+05   4.704314e-01   2.578171e+00      22344
       14    7.720375e+05   4.261541e+05   4.480137e-01   2.784407e+00      24130
       15    7.575438e+05   4.360300e+05   4.244162e-01   2.987403e+00      25916
       16    7.575438e+05   4.448043e+05   4.128336e-01   3.198021e+00      27702
       17    7.493192e+05   4.393746e+05   4.136350e-01   3.408740e+00      29488
       18    7.442237e+05   4.461845e+05   4.004699e-01   3.622386e+00      31274
       19    7.442237e+05   4.606246e+05   3.810671e-01   3.833010e+00      33060
       20    7.442237e+05   4.684432e+05   3.705613e-01   4.048867e+00      34846
       21    7.376148e+05   4.684432e+05   3.649216e-01   4.266168e+00      36632
       22    7.216422e+05   4.678996e+05   3.516182e-01   4.487643e+00      38418
       23    7.184630e+05   4.790798e+05   3.331879e-01   4.704879e+00      40204
       24    7.184630e+05   4.851569e+05   3.247295e-01   4.921035e+00      41990
       25    7.184630e+05   4.855915e+05   3.241245e-01   5.137653e+00      43776
       26    6.598698e+05   4.738143e+05   2.819580e-01   5.361330e+00      45562
       27    6.404863e+05   4.855915e+05   2.418393e-01   5.585230e+00      47348
       28    6.337599e+05   4.851738e+05   2.344517e-01   5.811326e+00      49134
       29    6.337599e+05   4.967139e+05   2.162428e-01   6.033616e+00      50920
       30    6.293701e+05   5.050416e+05   1.975444e-01   6.260361e+00      52706
       31    6.088077e+05   4.967139e+05   1.841203e-01   6.491571e+00      54492
       32    6.028812e+05   5.050416e+05   1.622868e-01   6.733449e+00      56278
       33    5.962405e+05   5.050416e+05   1.529566e-01   6.961139e+00      58064
       34    5.877859e+05   5.050416e+05   1.407729e-01   7.202975e+00      59850
       35    5.856300e+05   5.094496e+05   1.300828e-01   7.428922e+00      61636
       36    5.856300e+05   5.111242e+05   1.272233e-01   7.655393e+00      63422
       37    5.856300e+05   5.135643e+05   1.230568e-01   7.883653e+00      65208
       38    5.856300e+05   5.139124e+05   1.224624e-01   8.110602e+00      66994
       39    5.856300e+05   5.143420e+05   1.217288e-01   8.341516e+00      68780
       40    5.856300e+05   5.163903e+05   1.182312e-01   8.577273e+00      70566
       41    5.856300e+05   5.175452e+05   1.162592e-01   8.812001e+00      72352
       42    5.856300e+05   5.175475e+05   1.162552e-01   9.047505e+00      74138
       43    5.856300e+05   5.175475e+05   1.162552e-01   9.287532e+00      75924
       44    5.856300e+05   5.187268e+05   1.142414e-01   9.525509e+00      77710
       45    5.856300e+05   5.187268e+05   1.142414e-01   9.760768e+00      79496
       46    5.856300e+05   5.187268e+05   1.142414e-01   9.991682e+00      81282
       47    5.856300e+05   5.187268e+05   1.142414e-01   1.022682e+01      83068
       48    5.856300e+05   5.187268e+05   1.142414e-01   1.046331e+01      84854
       49    5.847353e+05   5.070967e+05   1.327756e-01   1.070000e+01      86640
       50    5.690759e+05   5.187268e+05   8.847517e-02   1.094639e+01      88426
       51    5.572637e+05   5.192932e+05   6.813729e-02   1.119262e+01      90212
       52    5.572637e+05   5.192932e+05   6.813729e-02   1.144323e+01      91998
       53    5.480985e+05   5.194214e+05   5.232101e-02   1.169302e+01      93784
       54    5.480985e+05   5.194896e+05   5.219670e-02   1.193647e+01      95570
       55    5.480985e+05   5.194896e+05   5.219670e-02   1.218329e+01      97356
       56    5.480985e+05   5.194896e+05   5.219670e-02   1.243950e+01      99142
       57    5.480985e+05   5.194896e+05   5.219670e-02   1.268879e+01     100928
       58    5.480985e+05   5.194896e+05   5.219670e-02   1.293691e+01     102714
       59    5.480985e+05   5.194896e+05   5.219670e-02   1.317849e+01     104500
       60    5.480985e+05   5.194896e+05   5.219670e-02   1.342800e+01     106286
       61    5.480985e+05   5.195993e+05   5.199641e-02   1.368159e+01     108072
       62    5.480985e+05   5.195993e+05   5.199641e-02   1.392741e+01     109858
       63    5.480985e+05   5.195993e+05   5.199641e-02   1.418535e+01     111644
       64    5.389053e+05   5.088605e+05   5.575157e-02   1.444954e+01     113430
       65    5.340461e+05   5.204619e+05   2.543631e-02   1.470719e+01     115216
       66    5.340461e+05   5.214920e+05   2.350758e-02   1.495451e+01     117002
       67    5.340461e+05   5.255546e+05   1.590024e-02   1.520475e+01     118788
       68    5.340461e+05   5.255546e+05   1.590024e-02   1.545330e+01     120574
       69    5.340461e+05   5.255546e+05   1.590024e-02   1.570904e+01     122360
       70    5.247515e+05   5.255546e+05  -1.530477e-03   1.597255e+01     124146
Converged. 70 iterations in 15.973 seconds.
