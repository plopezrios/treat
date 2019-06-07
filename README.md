TREAT
=====
Tail-regression estimator analysis toolbox.

This utility implements the statistical method described in
https://arxiv.org/abs/1903.07684 to analyze samples distributed according to
heavy-tailed distributions of known tail indices.

TREAT is a stand-alone Fortran utility using the MPI and LAPACK libraries.  The
code can be compiled with:

  mpif90 -o treat treat.sh -llapack

The compile.sh script is provided for reference, and also shows how to use
Intel's MKL LAPACK library.

TREAT can be run with, e.g.,

  ./treat
  mpirun -np 4 ./treat
  mpirun -np 20 ./treat < input > output &

TREAT uses a command-line interface which can be scripted easily, and the main
output entities are easily grep-able (notice the lines starting ASSESS, STAT,
and EVAL in the example below).  Type 'help' at the prompt for instructions.

Usage example
=============
The following should summarize the main capabilities of the code:


```
$ cat input
# Enable input echo.
set echo
# Generate sample of 10000 independent random numbers identically distributed
# according to heavy-tailed dstribution with principal exponent mu=2.5 and
# exponent increment dmu=0.5.
generate 10000 h(2.5) 0.5*h(3)
# Verify absence of serial correlation (independent random number assumption).
report corr
# Report basic statistical parameters (mean, variance, median, max, min).
report stats
# Perform tail-index estimation to verify mu.
assess mu at mlogq=1:7:13
# Set mu and dmu to known value.
set mu 2.5
set dmu 0.5
# Set constraint that leading-order coefficient is the same on both tails.
set constraint l1 = r1
# Automatically choose number of parameters and fit anchor point from
# provided ranges.
assess nparam,anchor at mlogq=1:7:13 for nparam=2:8 apply
# Evaluate tail-regression estimator of mean; while at it, make plots of the
# probability distribution (a Gaussian-kernel-smoothed histogram) and of the
# tail fits.
evaluate TRE plot fit to fit.plot plot P to P.plot
$ mpirun -np 4 ~/dev/treat/treat < input > out
$ cat out
==================================================
TREAT - tail-regression estimator analysis toolbox
==================================================

4 MPI processes available.

Type "help" for a list of commands.

TREAT> 
TREAT> 
Enabled input echo.

TREAT> 
TREAT> 
TREAT> 
TREAT> generate 10000 h(2.5) 0.5*h(3)

Generating sample drawn from model distribution
===============================================
Mean              :   0.000000000000E+00
Variance          : (undefined)
Asymptote         : 0.2523 * A^-2.5 + 0.1378 * A^-3
Number of samples : 10000

TREAT> 
TREAT> report corr

Serial correlation analysis
===========================
* Error factor for 1 chunks:
  ndecorr            1            2            4            8           16
   0-100%       0.9785       0.9879       0.9925       0.9850       1.0232

* Error factor for 3 chunks:
  ndecorr            1            2            4            8           16
   0- 33%       0.9611       0.9876       0.9270       1.0241       1.0187
  33- 67%       0.9171       1.2761       1.1302       1.1738       0.9585
  67-100%       0.9793       1.0277       0.9705       1.2786       1.0450

* Error factor for 5 chunks:
  ndecorr            1            2            4            8           16
   0- 20%       0.9409       0.9887       0.9292       0.9889       0.9618
  20- 40%       0.9164       0.9328       0.8482       0.9826       0.7426
  40- 60%       1.0677       1.0108       1.0675       1.1292       0.9483
  60- 80%       1.1376       1.0847       1.0870       1.1726       1.0715
  80-100%       0.9439       0.9247       1.0584       0.9029       0.9426

* Error factor for 10 chunks:
  ndecorr            1            2            4            8           16
   0- 10%       0.9136       0.8857       0.9503       0.9628       1.1435
  10- 20%       0.9795       0.8897       0.8802       0.8349       0.9736
  20- 30%       0.9026       0.8931       1.2579       1.0580       0.9589
  30- 40%       1.0251       0.9250       1.3103       0.7634       1.4027
  40- 50%       1.0072       0.9016       0.8953       0.7032       0.6515
  50- 60%       1.0099       1.0019       1.3154       0.8717       0.8264
  60- 70%       1.0944       1.2388       1.0366       0.7992       0.8001
  70- 80%       1.0638       1.1974       1.0023       0.9696       0.8866
  80- 90%       1.0513       0.8818       0.9839       0.7614       1.1924
  90-100%       0.8408       0.8769       1.3490       0.8381       1.9126

TREAT> 
TREAT> report stats

Basic statistics
================
STAT  Mean       -1.732806469400E-02   4.389419140337E-02
STAT  Variance    1.926700038956E+01   4.122410473299E+00
STAT  Median     -1.610182464356E-02
STAT  Max+        1.071345944405E+02
STAT  Max-       -1.634363045692E+02

TREAT> 
TREAT> assess mu at mlogq=1:7:13

Tail-index estimation
=====================
* Left tail:
  * mu    =    2.528347336752E+00   5.495068184818E-02 (=~ 2.50)
  * c     =    3.649785598566E-01   2.466284612679E-02
  * mlogq =    2.500609310290E+00   4.004632088884E-15
  * logA  =    6.229603133072E-01   2.406961053090E-02
* Right tail:
  * mu    =    2.503796303686E+00   5.566489847879E-02 (=~ 2.50)
  * c     =    3.566223337307E-01   2.416370164091E-02
  * mlogq =    2.499391060743E+00   2.669754725922E-15
  * logA  =    6.311245479142E-01   2.454157608187E-02

TREAT> 
TREAT> set mu 2.5

Set mu for both tails to 2.5.

TREAT> set dmu 0.5

Set dmu for both tails to 0.5.

TREAT> 
TREAT> set constraint l1 = r1

Defined constraint #1.

TREAT> 
TREAT> 
TREAT> assess nparam,anchor at mlogq=1:7:13 for nparam=2:8 apply

Fit assessment
==============
 -mlogq of left anchor-- -mlogq of right anchor- --1/A of left anchor--- --1/A of right anchor-- --nparam--- ------Norm estimator------- ------Mean estimator------- ----Variance estimator----- -----------Chi^2----------- RDC+0U12X
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00     0     0  1.000000E+00  0.000000E+00 -2.165397E-02  4.319924E-02  1.905599E+01  3.880868E+00  0.000000E+00  0.000000E+00 ---------
  1.0001E+00  8.8992E-16  1.0001E+00  8.8992E-16  2.9221E+00  1.0059E-01  2.9579E+00  1.2543E-01     6     6  9.997399E-01  3.207776E-04 -1.410710E-02  2.614265E-02  1.023490E-02  3.905993E-04  1.760221E-07  4.220414E-08 >>>PPPP-P
  1.4999E+00  1.1124E-15  1.4999E+00  1.1124E-15  1.2591E+00  2.7681E-02  1.2399E+00  2.6247E-02     4     4  9.999062E-01  4.391651E-04 -1.028418E-02  2.407995E-02  1.044007E-01  2.591598E-03  1.267127E-06  3.809136E-07 >>>PPPP-P
  1.9999E+00  2.4473E-15  1.9999E+00  2.4473E-15  7.8223E-01  1.4359E-02  7.7750E-01  1.4970E-02     4     4  9.996824E-01  3.745025E-04 -1.226231E-02  2.550397E-02  2.887123E-01  6.503963E-03  3.109104E-06  8.864851E-07 >>>PPPP-P
  2.5004E+00  1.7798E-15  2.5004E+00  1.7798E-15  5.3467E-01  1.3514E-02  5.3216E-01  1.2137E-02     4     4  9.998837E-01  3.275182E-04  1.814936E-02  3.060824E-02  5.456858E-01  1.254158E-02  6.201052E-06  1.587815E-06 >>>PPPP-P
  3.0007E+00  1.7798E-15  3.0007E+00  1.7798E-15  3.7518E-01  8.9397E-03  3.8004E-01  1.0713E-02     4     4  9.999019E-01  3.070510E-04 -3.511584E-03  3.247985E-02  8.654692E-01  2.155189E-02  1.111239E-05  3.413393E-06 >>>PPPP-P
  3.5016E+00  2.2248E-15  3.5016E+00  2.2248E-15  2.7759E-01  9.7757E-03  2.7351E-01  9.2454E-03     4     4  9.998950E-01  2.751463E-04 -2.106633E-02  3.587424E-02  1.236211E+00  3.279149E-02  1.717407E-05  5.054936E-06 >>>PPPP-P
  3.9981E+00  3.5597E-15  3.9981E+00  3.5597E-15  2.0684E-01  8.0787E-03  1.9299E-01  9.7953E-03     4     4  9.999301E-01  2.193870E-04 -3.028202E-02  4.041944E-02  1.666800E+00  5.125732E-02  3.234431E-05  1.215159E-05 >>>PPPP-P
  4.4963E+00  5.3395E-15  4.4963E+00  5.3395E-15  1.5005E-01  1.0183E-02  1.4346E-01  7.0853E-03     3     3  9.999711E-01  2.171437E-04 -1.939943E-02  3.153427E-02  2.157530E+00  7.961046E-02  7.073751E-05  2.638281E-05 >>>PPPP-P
  4.9982E+00  0.0000E+00  4.9982E+00  0.0000E+00  1.0819E-01  9.0184E-03  1.0834E-01  9.8599E-03     3     3  9.999702E-01  1.823568E-04 -1.413418E-02  3.604735E-02  2.704973E+00  1.216864E-01  1.266765E-04  5.227750E-05 >>>PPPP-P
  5.5090E+00  6.2294E-15  5.5090E+00  6.2294E-15  7.5536E-02  7.5321E-03  7.8040E-02  8.9613E-03     3     3  9.998841E-01  1.265195E-04 -1.067186E-02  4.539449E-02  3.367451E+00  1.958852E-01  2.251602E-04  1.057203E-04 >>>PPPP-P
  6.0117E+00  1.7798E-15  6.0117E+00  1.7798E-15  5.1237E-02  8.9803E-03  5.1686E-02  8.5871E-03     2     2  9.999714E-01  1.371379E-04 -1.296767E-02  3.736315E-02  4.220645E+00  3.527535E-01  1.021641E-03  7.882806E-04 >>>PPPP-P
  6.4695E+00  2.6698E-15  6.4695E+00  2.6698E-15  3.4649E-02  7.3265E-03  3.6118E-02  7.6665E-03     2     2  1.000009E+00  1.701706E-04 -1.327669E-02  4.466556E-02  5.316622E+00  6.192926E-01  1.869295E-03  1.130780E-03 >>>PPPP-P
  6.9590E+00  4.4496E-15  6.9590E+00  4.4496E-15  2.4736E-02  4.8830E-03  2.5061E-02  5.3386E-03     2     2  9.999433E-01  1.121600E-04 -1.087093E-02  4.267255E-02  6.887077E+00  9.780417E-01  2.464724E-03  1.692614E-03 >>>PPPP-P

Suggested anchor and nparam:
ASSESS  anchor  q   2.2313016014843E-01  q   2.2313016014843E-01
ASSESS  anchor  mlogq   1.5000000000000E+00  mlogq   1.5000000000000E+00
ASSESS  nparam     4   4
Setting anchor and nparam to suggested values.

TREAT> 
TREAT> 
TREAT> 
TREAT> evaluate TRE plot fit to fit.plot plot P to P.plot

Tail-regression estimation
==========================
* Standard estimator:
  * Mean                :  -1.828833062764E-02   4.256659091789E-02
  * Variance            :   1.947873528561E+01   3.782317526553E+00
* Tail-regression estimator:
  * Ac                  :  -1.610182464356E-02
  * Tail 1:
    * Fit anchor:
      * A               :   7.759258824663E-01   1.791018401722E-02
      * log|A-Ac|       :  -2.334133919930E-01   2.260236789134E-02
      * |A-Ac|^-1       :   1.263224760724E+00   2.854946916893E-02
      * M               :   2.232000000000E+03   0.000000000000E+00
      * q               :   2.232000000000E-01   0.000000000000E+00
      * -logq           :   1.499911087907E+00   1.112397802468E-15
    * Coefficients:
      * c_1             :   1.085447946474E-01   1.019343749880E-01
      * c_2             :   1.237847087835E+00   6.151848919350E-01
      * c_3             :  -1.787011769585E+00   1.024815616920E+00
      * c_4             :   5.497537363157E-01   5.082190568166E-01
    * Contrib. to norm  :   2.232215502867E-01   2.829831617539E-04
    * Contrib. to mean  :  -5.826036329011E-01   3.139614354294E-02
  * Tail 2:
    * Fit anchor:
      * A               :   7.909760521216E-01   1.463789679031E-02
      * log|A-Ac|       :  -2.144990013005E-01   1.814256324303E-02
      * |A-Ac|^-1       :   1.239444115296E+00   2.249946869875E-02
      * M               :   2.232000000000E+03   0.000000000000E+00
      * q               :   2.232000000000E-01   0.000000000000E+00
      * -logq           :   1.499911087907E+00   1.112397802468E-15
    * Coefficients:
      * c_1             :   1.085447946474E-01   1.019343749880E-01
      * c_2             :   1.377301688321E+00   6.176264601586E-01
      * c_3             :  -2.190220572644E+00   1.057902903265E+00
      * c_4             :   8.204950475582E-01   5.427306073092E-01
    * Contrib. to norm  :   2.230747553311E-01   2.877252624472E-04
    * Contrib. to mean  :   5.847078755369E-01   3.229229795009E-02
  * Centre:
    * Contrib. to norm  :   5.536000000000E-01   0.000000000000E+00
    * Contrib. to mean  :  -7.846823051187E-03   5.699129803671E-03
    * Contrib. to var.  :   1.042149756933E-01   2.558765908128E-03
  * Norm                :   9.998963056178E-01   3.611214539491E-04
  * Mean                :  -5.742580415377E-03   2.479621926337E-02

Plotted P to "P.plot".

Plotted fit to "fit.plot".

EVAL  Std.mean   -1.828833062764E-02   4.256659091789E-02
EVAL  Std.var.    1.947873528561E+01   3.782317526553E+00
EVAL  Norm        9.998963056178E-01   3.611214539491E-04
EVAL  Mean       -5.742580415377E-03   2.479621926337E-02

TREAT> 
Quitting.
```
