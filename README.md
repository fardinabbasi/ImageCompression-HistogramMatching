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
### SVD & FFT Method for Color Images
| Method | Comperession Rate = 0.1% | Comperession Rate = 0.5% | Comperession Rate = 1% | Comperession Rate = 4% | Comperession Rate = 8% | Comperession Rate = 10%| Comperession Rate = 12% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SVD | <img src="/readme_images/cs1.png"> | <img src="/readme_images/cs2.png"> | <img src="/readme_images/cs3.png"> | <img src="/readme_images/cs4.png"> | <img src="/readme_images/cs5.png"> | <img src="/readme_images/cs6.png"> | <img src="/readme_images/cs7.png"> |
| FFT | <img src="/readme_images/cf1.png"> | <img src="/readme_images/cf2.png"> | <img src="/readme_images/cf3.png"> | <img src="/readme_images/cf4.png"> | <img src="/readme_images/cf5.png"> | <img src="/readme_images/cf6.png"> | <img src="/readme_images/cf7.png"> |

In general, the FFT method provides **higher-quality** compression compared to the SVD method. This is due to the differences in their algorithms, where the FFT method allows for accurate **low-pass** and **high-pass** filtering. [Read More](https://ieeexplore.ieee.org/document/7424148)
## Denoising
Since noise in an image is like a **DC offset** in its **singular values**, it can be eliminated by calculating the differences between two singular values and cutting those whose differences are less than the **threshold**.
Here is the implementation of image denoising using the SVD method **from scratch**.
```ruby
def SVD_Denoise(filename, rank):
  img = Image.open(filename)
  img = np.asarray(img)
  denoised_img = np.zeros(img.shape)
  for rgb in range(img.shape[2]):
    U, S, V = np.linalg.svd(img[:,:,rgb])
    denoised_img[:,:,rgb] = np.matmul(np.matmul(U[:, :rank] , np.diag(S[:rank])) , V[:rank, :])

  for ind1, row in enumerate(denoised_img):
    for ind2, col in enumerate(row):
      for ind3, value in enumerate(col):
        if value < 0:
          denoised_img[ind1,ind2,ind3] = abs(value)
        if value > 255:
          denoised_img[ind1,ind2,ind3] = 255  
  return denoised_img.astype(np.uint8)
```
<img src="/readme_images/s_noise.png">

Because noises have **high frequencies**, the FFT of the image can be used to determine the **cut-off frequency** for preserving only a specific frequency range.
Here is the implementation of image denoising using the FFT method **from scratch**.

```ruby
def FFT_Denoise(filename, r):
  img = Image.open(filename)
  img = np.asarray(img)
  denoised_img = np.zeros(img.shape)

  for rgb in range(img.shape[2]):
    fft_img = np.fft.fft2(img[:,:,rgb])
    rows, cols = fft_img.shape
    fft_img[int(rows*r):int(rows*(1-r)),:] = 0
    fft_img[:, int(cols*r):int(cols*(1-r))] = 0
    denoised_img[:,:,rgb] = np.fft.ifft2(fft_img).real

  for ind1, row in enumerate(denoised_img):
    for ind2, col in enumerate(row):
      for ind3, value in enumerate(col):
        if value < 0:
          denoised_img[ind1,ind2,ind3] = abs(value)
        if value > 255:
          denoised_img[ind1,ind2,ind3] = 255  
  return denoised_img.astype(np.uint8)
```
<img src="/readme_images/f_noise.png">

## Histogram Matching
Histogram matching is a quick and easy way to "**calibrate**" one image to match another. In mathematical terms, it's the process of transforming one image so that the **cumulative distribution function** (CDF) of values in each band matches the CDF of bands in another image.

Here is the implementation of Histogram matching **from scratch** with the results.
```ruby
def Matching_Histogram(Reference, Source):
  num_bins = 255
  matched_image = Source
  for rgb in range(Source.shape[2]):
    # Calculate CDF of the images:
    ref_CDF, bins = CDF(Reference, rgb, num_bins)
    src_CDF, bins = CDF(Source, rgb, num_bins)
    # Normalizing the CDFs:
    ref_CDF = Normalize(ref_CDF)
    src_CDF = Normalize(src_CDF)
    # Matching the images:
    new_src = np.interp(Source[:,:,rgb].flatten(), bins[:-1], src_CDF)
    changed_src = np.interp(new_src, ref_CDF, bins[:-1])
    matched_image[:,:,rgb] = changed_src.reshape((Source.shape[0],Source.shape[1]))
  return matched_image
```
```ruby
def CDF(input, rgb, num_bins):
  Hist, bins = np.histogram(input[:,:,rgb].flatten(), num_bins)
  cdf = np.cumsum(Hist)
  return cdf, bins

def Normalize(input_CDF):
  return (255*input_CDF/input_CDF[-1]).astype(np.uint8)
```
<img src="/readme_images/matched.png">

## Modified Gram-Schmidt
To determine Q in **QR decomposition**, the Gram-Schmidt method is commonly employed. However, the traditional Gram-Schmidt method is **susceptible to rounding errors** and other issues. As an alternative, the **modified Gram-Schmidt** method can be utilized. The code provided below demonstrates the implementation of QR decomposition from scratch using both the Gram-Schmidt and modified Gram-Schmidt algorithms:
```ruby
def QR(A):
	r, c = A.shape
	Q = np.zeros((r, c),dtype=np.float64) # initialize matrix Q
	u = np.zeros((r, c),dtype=np.float64) # initialize matrix u
	u[:, 0] = copy.copy(A[:, 0])
	Q[:, 0] = u[:, 0] / np.linalg.norm(u[:, 0])
	for i in range(1, c):
		u[:, i] = A[:, i]
		for j in range(i):
			u[:, i] -= np.dot(A[:, i] , Q[:, j]) * Q[:, j] # get each u vector
		Q[:, i] = u[:, i] / np.linalg.norm(u[:, i]) # compute each e vetor
   # QT=np.transpose(Q)
	R = np.zeros((r, c),dtype=np.float64)
	for i in range(n):
		for j in range(i, c):
			R[i, j] = np.dot(A[:, j] , Q[:, i])
    #R=np.matmul(QT,A)
	return Q,R
```
```ruby
def QR_Modified_Decomposition(A):
	r, c = A.shape # get the shape of A
	Q = np.zeros((r, c),dtype=np.float64) # initialize matrix Q
   # u = np.zeros((n, m),dtype=np.float64) # initialize matrix u
	R = np.zeros((r, c),dtype=np.float64)
	u = copy.copy(A)
	for i in range(c):
		R[i,i]=np.linalg.norm(u[:, i])
		Q[:, i] = u[:, i] / np.linalg.norm(u[:, i])
		for j in range(i,n):
			R[i,j]= np.dot(Q[:, i] , u[:, j])
			u[:, j] -= (np.dot(Q[:, i] , u[:, j]))* Q[:, i] # get each u vector
	return Q,R
```
<img src="/readme_images/gram2.png">

According to the presented results, in larger row sizes, the modified Gram-Schmidt algorithm performs better and produces results similar to singular values of A. In contrast, classic Gram-Schmidt produces a higher error.
