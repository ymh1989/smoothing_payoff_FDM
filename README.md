##Smoothing Payoff for Finite Difference Method##

###Introduction###

It is difficult for derivatives with discontinuous payoff to be obtained a accurate numerical solution at around its maturity. This repo contains an implementation of methods of smoothing discontinous payoff for finite difference method(FDM). I consider the simulation for European digital call option which has a analytic solution\[1\]. 

###Methods###

- In general, it is known that the strike price should lie midway betweeen two grids points for increases the accuracy of the FDM[2]. The detailed descriptions about shifting the grid can be found in [3].

- Also, one of the method for increasing acccuracy of discontinous payoff on FDM is to smooth the payoff by using articial adjustment of payoff structure. This method is to supplement several grid points around strike price which has linear or nonlinear form payoff value.


###Environment###

- CPU : Intel(R) Core(TM) i5-6400 @ 2.7GHZ
- RAM : DDR3L 16GB PC3-12800
- Windows 10 Enterprise 64bit
- [Python 2.7](https://www.python.org/), [numpy 1.10.4](http://www.numpy.org/)
- [MATLAB 2016a](www.mathworks.com/products/matlab/)



###Numerical results###
- In this repo, I compare the RMSE(Root Mean Square Error) on interest area and maximum error, respectively, between analytic and numerical solution each three version.
- `version 0` : Discontious payoff (original method)
- `version 1` : Linear line
- `version 2` : Nonlinear curve 1
- `version 3` : Nonlinear curve 2

####1. Parameters####
- Test case : European digital call option
- FDM method : Fully implicit FDM
- You can modifiy the parameters.

|            | Stock | Strike | Maturity  | Riskless  <p>interest rate</p> | Volatility | Number of simulations | Cash  |
|------------|-------|--------|-----------|-------------------------|------------|-----------------------|-------|
| Parameters | 100.0 | 100.0  | 1.0/365.0 | 0.03                    | 0.3        | 10, 10, 10            | 100.0 |

####2. Results####
- RMSE is defined on 80.0~120.0.
- Maximum error is obtained on all of compuational domain.
- The figures are contained in `Code` directory.

|           | RMSE       | Maxerr     |
|-----------|------------|------------|
| version 0 | 0.02485242 | 0.07421257 |
| version 1 | 0.01629715 | 0.04338666 |
| version 2 | 0.01859031 | 0.04798860 |
| version 3 | 0.01535285 | 0.03539024 |


###Future work###
- In fact, this implementation is for increasing accuracy of ELS(Equity Linked Securities) which is a popular financial instrument in South Korea. I'm going to apply this method in pricing ELS for accurate pricing.
- For simplifying, I perform the experiments by using a fully implcit method. However, Crank-Nicolson(CN) method which incorporates both explicit and implicit method is a more accurate method since its second order accuracy and stability. If CN method could be applied in practical method, numerical errror will be decreased. Unfortunately, CN method exhibit undesirable qualities if the time step is relatively large to spatial step[2] and [4]. To treat the solution, Rannacher time stepping was introduced. Details of the algorithm can be found in references [5-6].


###Note###
- If you're interested in my works, please visit my [homepage](https://sites.google.com/site/yoomh1989/).

###Reference###

\[1\] Haug, Espen Gaarder. The complete guide to option pricing formulas. McGraw-Hill Companies, 2007.
\[2\] Tavella, Domingo, and Curt Randall. Pricing financial instruments: The finite difference method. Vol. 13. John Wiley & Sons, 2000.
\[3\] Pooley, David M., Kenneth R. Vetzal, and Peter A. Forsyth. "Convergence remedies for non-smooth payoffs in option pricing." Journal of Computational Finance 6.4 (2003): 25-40.
\[4\] Jeong, Darae, and Junseok Kim. "A comparison study of ADI and operator splitting methods on option pricing models." Journal of Computational and Applied Mathematics 247 (2013): 162-171.
\[5\] Rannacher, Rolf. "Finite element solution of diffusion problems with irregular data." Numerische Mathematik 43.2 (1984): 309-327.
\[6\] Giles, Michael B., and Rebecca Carter. "Convergence analysis of Crank-Nicolson and Rannacher time-marching." (2005).