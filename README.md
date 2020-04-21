# RDDP_PowerSystem
Robust Dual Dynamic Programming for solving real-time operation problems in power system.

The RDDP solves multi-stage minimax problems in the following form:

<img src="/tex/16482209de72951a8693e6f3a40bce97.svg?invert_in_darkmode&sanitize=true" align=middle width=197.46581579999997pt height=27.91243950000002pt/>
subject to <img src="/tex/c7a6e236b62ff0e22574281733ee86f9.svg?invert_in_darkmode&sanitize=true" align=middle width=77.61033389999999pt height=22.831056599999986pt/>
<p align="center"><img src="/tex/949902a141187713e75fd0734c80a936.svg?invert_in_darkmode&sanitize=true" align=middle width=62.45248844999999pt height=14.2073316pt/></p>

In which <img src="/tex/6f7c377862afd23cbee482bbb9658658.svg?invert_in_darkmode&sanitize=true" align=middle width=68.18504714999999pt height=24.65753399999998pt/> is called the stage <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> worst-case cost-to-go function, defined as below:

<img src="/tex/6b654cdf3e85844babac2ffd08002644.svg?invert_in_darkmode&sanitize=true" align=middle width=291.83060279999995pt height=24.65753399999998pt/>

<img src="/tex/85076c56634338e869c856fdcaf963b2.svg?invert_in_darkmode&sanitize=true" align=middle width=550.81504335pt height=93.84671339999997pt/>

The worst-case cost-to-go function can be bounded from above and from below.

1. **Upper Bound**

<img src="/tex/f62fe982b84615153da9a60daf2acf86.svg?invert_in_darkmode&sanitize=true" align=middle width=66.31741544999998pt height=26.97711060000001pt/>
<img src="/tex/51b0862b4295783e6677244ec61af7a3.svg?invert_in_darkmode&sanitize=true" align=middle width=442.83160515pt height=67.54534379999997pt/>
subject to <img src="/tex/1dffb4e7b499ba4ecd6650b33eed1e2b.svg?invert_in_darkmode&sanitize=true" align=middle width=48.99537059999999pt height=22.831056599999986pt/>

2. **Lower Bound**

<img src="/tex/3ecc99fdcf1dffe238e04d68c7fb474c.svg?invert_in_darkmode&sanitize=true" align=middle width=463.2241119pt height=70.30311089999996pt/>


