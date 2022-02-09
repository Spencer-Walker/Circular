# File: mymodule.py
import mpmath as mp
import numpy as np

def yukawa(*args,**kwargs):
  n = args[0]
  i = args[1]
  l = args[2]
  k = args[3]
  Zc = args[4]
  g = args[5] 
  ret = complex( -k*np.sqrt(mp.gamma(n+1)*mp.gamma(i+1)/((n+l+1)*(i+l+1)*mp.gamma(n+2*l+2)*mp.gamma(i+2*l+2)))
    *Zc*(mp.gamma(n+i+2*l+2)/(mp.gamma(n+1)*mp.gamma(i+1)))*(g**(2*l+2)/(g+1)**(n+i+2*l+2))
    *mp.hyp2f1(-i,-n,-n-i-2*l-1,1-g**2))
  return ret


def shell(*args,**kwargs):
  n = args[0]
  i = args[1]
  l = args[2]
  g = args[3]
  a = args[4]
  if n != 0: 
    ret = complex(-np.sqrt(mp.gamma(n+1)*mp.gamma(i+1)/((n+l+1)*(i+l+1)*mp.gamma(n+2*l+2)*mp.gamma(i+2*l+2)))
      *(0.5*a)*( 2*(n+l+1)* (mp.gamma(n+i+2*l+2)/(mp.gamma(n+1)*mp.gamma(i+1)))
      *(g**(2*l+2)/(g+1)**(n+i+2*l+2))*mp.hyp2f1(-i,-n,-n-i-2*l-1,1-g**2)
      -(n+1)*(mp.gamma(n+i+2*l+3)/(mp.gamma(n+2)*mp.gamma(i+1)))*(g**(2*l+2)/(g+1)**(n+i+2*l+3))
      *mp.hyp2f1(-i,-n-1,-n-i-2*l-2,1-g**2) -(n+2*l+1)*(mp.gamma(n+i+2*l+1)/(mp.gamma(n)*mp.gamma(i+1)))
      *(g**(2*l+2)/(g+1)**(n+i+2*l+1))*mp.hyp2f1(-i,-n+1,-n-i-2*l,1-g**2)))
  else:
    ret = complex(-np.sqrt(mp.gamma(n+1)*mp.gamma(i+1)/((n+l+1)*(i+l+1)*mp.gamma(n+2*l+2)*mp.gamma(i+2*l+2)))
      *(0.5*a)*( 2*(n+l+1)* (mp.gamma(n+i+2*l+2)/(mp.gamma(n+1)*mp.gamma(i+1)))
      *(g**(2*l+2)/(g+1)**(n+i+2*l+2))*mp.hyp2f1(-i,-n,-n-i-2*l-1,1-g**2)
      -(n+1)*(mp.gamma(n+i+2*l+3)/(mp.gamma(n+2)*mp.gamma(i+1)))*(g**(2*l+2)/(g+1)**(n+i+2*l+3))
      *mp.hyp2f1(-i,-n-1,-n-i-2*l-2,1-g**2)))
  return ret



