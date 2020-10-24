# golem-parallel-matplotlib

This is my submission for a [gitcoin bounty](https://gitcoin.co/issue/golemfactory/yagna/703/100023964) regarding the use of Golem for concurrent multiprocessing.

For this project, various statistical analyses are performed on circadian rhythm measurements in human test subjects.  The anonymized test data can be found in [datasets/ppd](datasets/ppd).  A nonlinear least squares regression is performed to fit the data to a cosine curve, and various calculations are performed with respect to the MESOR (Midline Estimating Statistics Of Rhythm).  The analysis is parameterized by a mesor threshold percentage.  For this demonstration, each dataset was run 3 times, with the thresholds set at 25%, 50%, and 75% respectively.  The analysis algorithm is under active development; we noticed that currently 3 datasets cause exceptions at `threshold=25`, but since there is no way to abort failing tasks in Golem, we excluded them from the demo (see [failures](datasets/ppd/failures)).

YouTube video demo [here](https://youtu.be/hflrBq2OXwA).<br>
An animated .gif quickly cycling through the results is [here](received/mogged/animated.gif).
