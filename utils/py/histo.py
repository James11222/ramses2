#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse
import ramses as ram

parser = argparse.ArgumentParser()
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="specify a path")
parser.add_argument("--log", help="plot log variable",action="store_true")
parser.add_argument("--bkp", help="plot log variable",action="store_true")
parser.add_argument("--out", help="output a png image")
parser.add_argument("--prefix", help="specify a file prefix")
parser.add_argument("--min", help="specify a minimum variable value")
parser.add_argument("--max", help="specify a maximum variable value")
parser.add_argument("--var", help="specify a variable number")
parser.add_argument("--xcen", help="specify the image center x-coordinate")
parser.add_argument("--ycen", help="specify the image center y-coordinate")
parser.add_argument("--zcen", help="specify the image center z-coordinate")
parser.add_argument("--rad", help="specify the image radius")
args = parser.parse_args()
# path the the file
path = args.path
prefix = args.prefix
ivar = args.var
vmin = args.min
vmax = args.max
radius = args.rad
xcenter = args.xcen
ycenter = args.ycen
zcenter = args.zcen
log = args.log
backup = args.bkp

if xcenter==None:
    xcenter=None
else:
    xcenter=float(xcenter)
if ycenter==None:
    ycenter=None
else:
    ycenter=float(ycenter)
if zcenter==None:
    zcenter=None
else:
    zcenter=float(zcenter)
if radius==None:
    radius=None
else:
    radius=float(radius)
center=np.array([xcenter,ycenter,zcenter])
log0=None
if log:
    log0=1
if ivar==None:
    ivar=0
else:
    ivar=int(ivar)-1
if prefix==None:
    prefix="hydro"
if path==None:
    path="./"
else:
    path=path+"/"
if prefix=="hydro":
    isort=0
if prefix=="peak_grid":
    isort=1
if prefix=="grav":
    isort=0
if prefix=="rt":
    isort=0

nout = args.nout
print("Reading output number ",nout)

c=ram.rd_cell(nout,path=path,prefix=prefix,center=center,radius=radius,backup=backup)
i=ram.rd_info(nout,path=path,backup=backup)
m_p=1.66e-24
k_b=1.38e-16
d=c.u[0]
if(backup):
    p=(i.gamma-1)*(c.u[4]-0.5*(c.u[1]**2+c.u[2]**2+c.u[3]**2)/c.u[0])
else:
    p=c.u[4]

plt.hist2d(np.log10(d*i.unit_d/m_p),np.log10(p/d/k_b*m_p*(i.unit_l/i.unit_t)**2),weights=c.u[0]*c.dx**3,density=True,bins=100,norm=colors.LogNorm())
plt.xlabel('log10(n_H) [H/cc]')
plt.ylabel('log10(T/mu) [K]')
plt.title('Mass-weighted Histogram with Logarithmic Color Scale')

if args.out:
    plt.savefig(args.out)

plt.show()
