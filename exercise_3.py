#######################################################################################
# 
#
#Exercise 3  of Homework-Part1
#
#Calculates and plots the Fourier transform of various functions.
#It was also made a comparion with theoretical values to see the quality of the aproximations
#introduced by each method.
#
#########################################################################################

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import math
import cmath

import scipy.signal
from scipy.fftpack import fft, fftshift, fftfreq, ifftshift, ifft


def build_f(a,times):
	
	vals=np.zeros(len(times))

	for i in range(len(times)):
		if abs(times[i])<a:
			vals[i]=1
	return vals


#alinea a 

def Ftransform_analit(t_max,a,step):
	#builds the f function FF calculated by numpy
	w_g=0

	freqs=np.arange(-t_max,t_max,step)

	#builds the analitical FF for f(t)
	vals=2*a*np.sinc(freqs*a)

	#builds the analitical FF for g(t)
	time_g=np.array([w])
	val_g=np.array([2*np.pi])

	#builds the box function for f 
	f_t,f=build_f(a,t_max,step)


	w_fg=1
	a_fg=1
	time_fg=np.array([w])
	val_fg=np.array([2*a_fg])




# FFT of f function

def Fourier_trans_f(a=5,size=100):

	
	N=2**12

	step=2*size/N

	#creates the time array for building the f function
	time=np.arange(-size,size,step)
	n=len(time)

	#creates the f function
	valor_f=build_f(a,time)

	#aplica um fftshift para centrar a frequencia zero aos valores calculados para a Fourier trasnform. M
	#multiplicado pelo step de forma a normalizar os valores

	fft=np.fft.fftshift(np.fft.fft(valor_f,n=N))*step


	#feito o mesmo shift ao array de frequencias para N pontos espaçados de step
	fft_time=np.fft.fftshift(np.fft.fftfreq(N,d=step))


	#devia ser 2pi f mas no sinc ja tem o pi
	#fftfreq retorna os valores em Hz, mas para a expressao analitica queremos uma freq angular
	w=2*fft_time
	certo_f=2*a*np.sinc(w*a)

	return fft


# FFT of g function

def Fourier_trans_g(w,size=100):
	N=2**12
	
	
	step=2*size/N

	#creates the time array for building the f function
	time=np.arange(-size,size,step)
	n=len(time)

	#creates the f function
	valor_g=np.exp(1j*w*time)


	fft=np.fft.fftshift(np.fft.fft(valor_g,n=N))*step

	
	fft_time=np.fft.fftshift(np.fft.fftfreq(N,d=step))

	
	time_g=np.array([w])
	val_g=np.array([2*np.pi])

	return fft

# FFT of f*g function

def Fourier_trans_fg(a=2,w=1,size=100):

	N=2**12
	
	
	step=2*size/N

	#creates the time array for building the f function
	time=np.arange(-size,size,step)
	n=len(time)

	#creates the functions
	valor_g=np.exp(1j*w*time)
	valor_f=build_f(a,time)
	
	fg=valor_f*valor_g


	#calculates the fft
	fft=np.fft.fftshift(np.fft.fft(fg,n=N))*step
	fft_time=np.fft.fftshift(np.fft.fftfreq(N,d=step))


	return fft


# FFT of the convolution of f with g

def convolution(a=5,w=1,size=100):
	N=2**12

	step=2*size/N

	#creates the time array for building the f function
	time=np.arange(-size,size,step)
	n=len(time)

	#creates the functions
	valor_g=np.exp(1j*w*time)
	valor_f=build_f(a,time)

	# calculates the fft for both functions
	ft_f=np.fft.fft(valor_g,n=N)
	ft_g=np.fft.fft(valor_f,n=N)


	val=np.multiply(ft_f,ft_g)
	
	#calculates the inverse FFT of the multiplication of the FFT of f with the FFT of g

	convol=np.fft.fftshift(np.fft.ifft(val))
	

	return convol



################################################################
'''

Calculates and plots the results for different values of w and a.


'''
################################################################

# values of a and W
plot=1
if plot:
	values_w=[1,4]
	values_a=[2,3]
else:
	values_w=[]
	values_a=[]

size=60
N=2**12
step=2*size/N

time=np.arange(-size,size,step)
fft_time=2*np.pi*np.fft.fftshift(np.fft.fftfreq(N,d=step))


for a in values_a:
	
	f=Fourier_trans_f(a)

	pl.figure('FFT da função f')
	pl.title('FFT da função f')
	pl.legend()

	
	pl.title('Amplitude')
	pl.plot(fft_time,np.abs(f),label='a= %s'%a)

	plt.gca().set_xlim(-size,size)
	pl.legend()


	pl.subplot(2,1,2)
	pl.title('Fase')
	pl.plot(fft_time,np.angle(f),label='a= %s'%a)
	plt.gca().set_xlim(-size,size)
	plt.tight_layout(rect=[0, 0.03, 1, 0.95])
	pl.legend()
	

	



	for w in values_w:
		fg=Fourier_trans_fg(a,w)
		conv=convolution(a,w)

		pl.figure('F*G')
		pl.title('Amplitude da transformada')
		pl.plot(fft_time,np.abs(fg),label='a= %s'%a+';w='+str(w))
		plt.gca().set_xlim(-size,size)
		pl.legend()

		pl.subplot(2,1,2)
		pl.plot(fft_time,np.angle(fg),label='a='+str(a)+';w='+str(w))
		plt.gca().set_xlim(-size,size)
		pl.legend()
	

		pl.figure('Convolution')
		#pl.subplot(2,1,1)
		pl.plot(time,np.abs(conv),label='a='+str(a)+';w='+str(w))
		plt.gca().set_xlim(-size,size)
		pl.legend()
		
		pl.subplot(2,1,2)
		pl.plot(time,np.angle(conv),label='a='+str(a)+';w='+str(w))
		plt.gca().set_xlim(-size,size)
		pl.legend()



for w in values_w:
	g=Fourier_trans_g(w)

	pl.figure('FFT da função g')
	pl.title('Amplitude da transformada')
	pl.plot(fft_time,np.abs(g),label='w=%s'%w)
	plt.gca().set_xlim(-5,5)

	pl.legend()
	pl.subplot(2,1,2)
	pl.plot(fft_time,np.angle(g),label='a='+str(a)+';w='+str(w))
	plt.gca().set_xlim(-size,size)
	pl.legend()




pl.show()