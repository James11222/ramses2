#!/usr/bin/env python3
import pyvista as pv
import ramses as ram
import numpy as np
import argparse

# Check if we should use non-interactive backend
import sys
if "--no-display" in sys.argv:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    print("Using non-interactive backend for batch processing")

parser = argparse.ArgumentParser()
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="specify a path")
parser.add_argument("--log", help="plot log variable",action="store_true")
parser.add_argument("--out", help="output a png image")
parser.add_argument("--prefix", help="specify a file prefix")
parser.add_argument("--col", help="choose the color map")
parser.add_argument("--min", help="specify a minimum variable value for colorbar")
parser.add_argument("--max", help="specify a maximum variable value for colorbar")
parser.add_argument("--var", help="specify a variable number")
parser.add_argument("--xcen", help="specify the image center x-coordinate")
parser.add_argument("--ycen", help="specify the image center y-coordinate")
parser.add_argument("--zcen", help="specify the image center z-coordinate")
parser.add_argument("--rad", help="specify the image radius")
parser.add_argument("--no-display", help="prevent GUI display (useful for batch processing)",action="store_true")

args = parser.parse_args()

path = args.path
prefix = args.prefix
ivar = args.var
vmin = args.min
vmax = args.max
col = args.col
radius = args.rad
xcenter = args.xcen
ycenter = args.ycen
zcenter = args.zcen
log = args.log
no_display = args.no_display

# Convert vmin and vmax to float if provided
if vmin is not None:
    vmin = float(vmin)
if vmax is not None:
    vmax = float(vmax)

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

log0=False
if log:
    log0=True

if ivar==None:
    ivar=1
else:
    ivar=int(ivar)

if prefix==None:
    prefix="hydro"
if path==None:
    path="./"
else:
    path=path+"/"

if col==None:
    col="jet"

nout = args.nout
print("Reading output number ",nout)

c=ram.rd_cell(nout,path=path,prefix=prefix,center=center,radius=radius)
x=c.x[0]
y=c.x[1]
z=c.x[2]
print(ivar,c.nvar)
if ivar <= c.nvar:
    xx = c.u[ivar-1]
elif ivar == 15: # temperature
    xx = c.u[4]/c.u[0]
elif ivar == 16: # magnetic energy
    xx = 0.5*(c.u[5]**2+c.u[6]**2+c.u[7]**2)
elif ivar == 17: # kinetic energy
    xx = 0.5*(c.u[1]**2+c.u[2]**2+c.u[3]**2)
else:
    print("unknown variable: use rho instead")
    xx = c.u[0]

print("min=",np.min(xx)," max=",np.max(xx))
min_val = 1e-3*np.max(xx)
data = ram.mk_cube(x,y,z,c.dx,xx)
grid = pv.ImageData()
grid.dimensions = np.array(data.shape) + 1
if log0:
    grid.cell_data["values"] = np.log10(data.flatten(order="F")+min_val)
else:
    grid.cell_data["values"] = data.flatten(order="F")

# Color range: --min/--max are given in linear units like in amr2img.py
values = grid.cell_data["values"]
if vmin is None:
    vmin = np.min(values)
elif log0:
    vmin = np.log10(vmin)
if vmax is None:
    vmax = np.max(values)
elif log0:
    vmax = np.log10(vmax)

pl = pv.Plotter(off_screen=no_display, window_size=[1600, 1600])
pl.add_volume(grid, scalars="values", cmap=col, opacity="sigmoid", clim=[vmin, vmax])
pl.add_bounding_box()

if args.out:
    pl.show(screenshot=args.out)
else:
    pl.show()
