# Efficient computation of first passage times in Kou's jump-diffusion model (Belkaid & Utzet, 2016)

This work aims to improve the performance of the algorithms developed by Belkaid and Utzet by translating them to a pure C numerical execution from symbolic Maple scripts.

## Speedup
Using `Digits:=30;` in Maple and 80-bit floats (`long double`) and 64-bit ints (`int64_t`) in C.
 
 - 1st function (f1+euler+suma): x723  (from 868ms to 1.20ms)
 - 2nd function (f2+euler+suma): x899  (from 1259ms to 1.40ms)

Using `Digits:=500;` in Maple and GNU GMP, MPFR and MPC multiprecision libraries in C.

 - 3rd function (f1+gaver+fn): x14.2 (from 3364ms to 236.94ms)

Tested in an Intel(R) Core(TM)2 Duo CPU E8400 @ 3.00GHz with 6GB of ram and Linux OS (Xubuntu 16.10 64bit).

## Order of accuracy
Computed comparing results for *n*=12 and *n*=50 for the two first functions:

- 1st function: result: 0.2558436487, order of accuracy: 10^(-8)
- 2nd function: result: 0.2236160448, order of accuracy: 10^(-7)
- 3rd function: results have all digits exactly computed.

## Stability of Gaver with respect to *n*
Stability with respect to *n* and time of computation for Gaver algorithm: see <a href="https://github.com/oscarsaleta/KouJDM/blob/master/results.md">results.md</a>.
