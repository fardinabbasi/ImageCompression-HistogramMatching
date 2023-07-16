# ImageCompression-HistogramMatching
## Lagrange Interpolation
**Manually** calculate the Lagrange interpolation for the given points.

$x = [0,0.6,1.03,1.39,1.76,2.09,2.29]$

$f(x) = [2.718,0.797,0.368,0.597,1.712,2.718,2.7]$

As an example, the calculation for the **quartic** function is provided below.

<img src="/readme_images/p4.jpg">

Implementing Lagrange interpolation **from scratch** for different degrees using the following code.
```ruby
def Lagrange (Lx, Ly,x):
	est_y=[];
	for i in x:
		y=0
		for k in range ( len(Lx) ):
			t=1
			for j in range ( len(Lx) ):
				if j != k:
					t=t* ((i-Lx[j]) /(Lx[k]-Lx[j]))
			y+= t*Ly[k]
		est_y.append(y)
	return est_y
```
The results are displayed below:

| Lagrange interpolation | Error |
| --- | --- |
| <img src="/readme_images/lagr_result.png"> | <img src="/readme_images/error.png"> |

The **best estimation** for n+1 points of data is n degree polynomial and in this case, the 6th degree has the lowest error.
## LU Decomposition

<h2> Part 1: Lagrangian Interpolation & LU Decomposition </h2>
<h2> Part 2: Using FFT & SVD for Image Denoising </h2>
<h2> Part 3: Histogram Matching </h2>
<h2> Part 4: Modified Gram-Schmidt </h2>
