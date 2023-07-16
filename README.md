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
Implementing LU Decomposition **from scratch** for an **upper Hessenberg matrix** using the provided code.
```ruby
def lu(A):
	[r,c] = np.shape(A)
	U = A.astype('float32')
	L = np.eye(r)
	for i in range(r-1):
		fac=U[i+1,i]/U[i,i]
		U[i+1,:]-=fac*U[i]
		L[i+1,i]=fac
	return L, U
```
$$A =
\begin{matrix}
1 & 4 & 2 & 3 \\
3 & 4 & 1 & 7 \\
0 & 2 & 3 & 4 \\
0 & 0 & 1 & 3 \\
\end{matrix}
\rightarrow lu(A): L = \begin{matrix}
1 & 0 & 0 & 0 \\
3 & 1 & 0 & 0 \\
0 & -0.25 & 1 & 0 \\
0 & 0 & 0.5714 & 1 \\
\end{matrix} ,  U = \begin{matrix}
1 & 4 & 2 & 3 \\
0 & -8 & -5 & -2 \\
0 & 0 & 1.75 & 3.5 \\
0 & 0 & 0 & 1 \\
\end{matrix}$$
## Image Compression
### SVD & FFT Method for Grayscale Images
| Method | Comperession Rate = 0.1% | Comperession Rate = 0.5% | Comperession Rate = 1% | Comperession Rate = 4% | Comperession Rate = 8% | Comperession Rate = 10%| Comperession Rate = 12% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SVD | <img src="/readme_images/s1.png"> | <img src="/readme_images/s2.png"> | <img src="/readme_images/s3.png"> | <img src="/readme_images/s4.png"> | <img src="/readme_images/s5.png"> | <img src="/readme_images/s6.png"> | <img src="/readme_images/s7.png"> |
| FFT | <img src="/readme_images/f1.png"> | <img src="/readme_images/f2.png"> | <img src="/readme_images/f3.png"> | <img src="/readme_images/f4.png"> | <img src="/readme_images/f5.png"> | <img src="/readme_images/f6.png"> | <img src="/readme_images/f7.png"> |

<h2> Part 2: Using FFT & SVD for Image Denoising </h2>
<h2> Part 3: Histogram Matching </h2>
<h2> Part 4: Modified Gram-Schmidt </h2>
