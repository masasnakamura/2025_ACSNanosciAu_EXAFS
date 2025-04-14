import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime
import larch
import shutil
from larch import Interpreter

session = Interpreter()

import tools

################
### Settings ###
################

### input start ###
sname = "GaPGM" ### Sample name
temps = ["008K", "050K", "100K", "150K", "200K", "300K", "100C"] ### Measurement temperature in str
tempsdict = {"008K": 8.4, "050K": 50.0, "100K": 100.0, "150K": 150.0, "200K": 200.0, "300K": 300.0, "100C": 373.0} ### Temperatures in float
tempcolors = {"008K": (193/255,0/255,253/255),"050K": (0/255,95/255,245/255), "100K": (0/255,183/255,180/255),"150K": (0/255,160/255,82/255),"200K": (183/255,143/255,0/255),"300K": (252/255,0/255,60/255),"100C":(255/255,0/255,177/255)} ### color settings for figure export
fitspace = "q" ### k, r, q, w(= wavelet)

elemsdict = {"PGM": ["Ru", "Rh", "Pd", "Ir", "Pt"], "GaPGM": ["Ru", "Rh", "Pd", "Ir", "Pt", "Ga"], "InPGM": ["Ru", "Rh", "Pd", "Ir", "Pt", "In"],"SnPGM": ["Ru", "Rh", "Pd", "Ir", "Pt", "Sn"]}
elemcolors = {"Ru": (252/255,0/255,60/255), "Rh": (0/255,95/255,245/255), "Pd": (183/255,143/255,0/255), "Ir": (0/255,183/255,180/255), "Pt": (193/255,0/255,253/255), "Ga": (234/255,115/255,0/255), "In": (126/255,161/255,0/255), "Sn": (125/255,121/255,127/255)} ### color settings for figure export
elems = elemsdict[sname] ### Absorption elements
elems_scat = elems ### Scattering elements
xtol = 4 ### Tolerance factor
ftol = 6 ### Tolerance factor
gtol = 6 ### Tolerance factor
export_dir = f"./out/{sname:s}" ### Export directory
### input end ###

### fitting conditions prepared in csv file ###
pars = pd.read_csv(f"./input_pars_{sname:s}.csv", header=0, index_col=0,
                   dtype={"edge": str, "composition": float, "amp": float, "enot_init": float, "r_init": float, "pre_min": float, "pre_max": float, "post_min": float, "post_max": float, "bkg_max": float, "kweight": int, "kmin_ftf": float, "kmax_ftf": float, "dk_ftf": float, "rmin_ftr": float, "rmax_ftr": float, "dr_ftr": float})

for elem1, elem2 in zip(sorted(elems), sorted(list(pars.index))):
    if elem1 != elem2:
        print("Check elements!!!")
        quit()

### prepare data in athena .prj file ###
dnames = np.empty((len(elems), len(temps)), dtype=object)
for i,elem in enumerate(elems):
    for j,temp in enumerate(temps):
        edge = pars["edge"][elem]
        dnames[i,j] = f"{sname:s}_{temp:s}_{elem:s}"

os.makedirs(export_dir, exist_ok=True)
os.makedirs(f"{export_dir:s}/dat", exist_ok=True)
os.makedirs(f"{export_dir:s}/fit", exist_ok=True)

outpath = f"{export_dir:s}"
shutil.copy(__file__, outpath)
shutil.copy(f"./input_pars_{sname:s}.csv", outpath)

#######################
### Data Processing ###
#######################

dats = np.empty_like(dnames, dtype=object)
normfs = {}

for i in range(dnames.shape[0]):
    file = larch.io.read_athena(f'./in/{sname:s}/{sname:s}_{elems[i]:s}{pars["edge"][elems[i]]:s}.prj')
    e0 = None
    for j in range(dnames.shape[1]):
        dname = dnames[i,j]
        elem = elems[i]
        dat = getattr(file, dname)

        tools.norm_frac_edge(dat, pars.loc[elem], e0=e0, elem=elem, edge=pars["edge"][elem])
        if j==0:
            e0 = dat.e0

        larch.xafs.autobk(dat, rbkg=pars["rbkg"][elem], e0=dat.e0, edge_step=dat.edge_step, kmin=0, kmax=pars["bkg_max"][elem], kweight=2, nclamp=2, clamp_lo=0, clamp_hi=0, group=dat)

        larch.xafs.xftf(dat, kweight=pars["kweight"][elem], kmin=pars["kmin_ftf"][elem], kmax=pars["kmax_ftf"][elem], dk=pars["dk_ftf"][elem], window="hanning", group=dat)

        larch.xafs.xftr(dat, rmin=pars["rmin_ftr"][elem], rmax=pars["rmax_ftr"][elem], dr=pars["dr_ftr"][elem], window="hanning", group=dat)

        outpath = f"{export_dir:s}/dat/{dname:s}.png"
        tools.fig_norm(dat, pars.loc[elem], save_flag=True, export_path=outpath)

        dats[i,j] = dat

    tools.fig_temps(export_dir, sname, elem, temps, tempcolors, save_flag=True)

###############################
### Settings before Fitting ###
###############################

dsets = np.empty_like(dats, dtype=object)
for i,abs in enumerate(elems):
    edge = pars["edge"][abs]
    tran = larch.xafs.feffit_transform(fitspace=fitspace, kmin=pars["kmin_ftf"][abs], kmax=pars["kmax_ftf"][abs], window="hanning", kw=pars["kweight"][abs], dk=pars["dk_ftf"][abs], rmin=pars["rmin_ftr"][abs], rmax=pars["rmax_ftr"][abs], dr=pars["dr_ftr"][abs], rwindow="hanning")
    for j,temp in enumerate(temps):
        fitpaths = []
        for scat in elems_scat:
            path = larch.xafs.feffpath(f"./feffinp/{abs:s}{edge:s}-{scat:s}/feff0001.dat")
            fitpaths.append(larch.xafs.feffpath(path.filename, degen=1, s02 = f"amp_{abs:s}*n_{abs:s}{scat:s}/x_{abs:s}", e0 = f"enot_{abs:s}", deltar = f"r_{abs:s}_{temp:s}+r_{scat:s}_{temp:s}-reff", sigma2 = f"ss_{abs:s}_{temp:s}+ss_{scat:s}_{temp:s}"))
        dsets[i,j] = larch.xafs.feffit_dataset(data=dats[i,j], pathlist=fitpaths, transform=tran)

exafspars = larch.fitting.param_group()
for i,elem in enumerate(elems):
    edge = pars["edge"][elem]
    exafspars.__setattr__(f"x_{elem:s}", larch.fitting.param(pars["composition"][elem]))
    exafspars.__setattr__(f"amp_{elem:s}", larch.fitting.param(pars["amp"][elem]))
    exafspars.__setattr__(f"enot_{elem:s}", larch.fitting.guess(pars["enot_init"][elem], vary=True))
    for j,temp in enumerate(temps):
        exafspars.__setattr__(f"r_{elem:s}_{temp:s}", larch.fitting.guess(pars["r_init"][elem], vary=True))
        exafspars.__setattr__(f"ss_{elem:s}_{temp:s}", larch.fitting.guess(larch.xafs.sigma2_eins(tempsdict[temp], 200., larch.xafs.feffpath(f"./feffinp/{elem:s}{edge:s}-{elem:s}/feff0001.dat"))/2., vary=True))
    for k,elem_scat in enumerate(elems_scat):
        if i<=k:
            exafspars.__setattr__(f"n_{elem:s}{elem_scat:s}", larch.fitting.guess(10.*pars["composition"][elem]*pars["composition"][elem_scat], min=0., vary=True))
        else:
            exafspars.__setattr__(f"n_{elem:s}{elem_scat:s}", larch.fitting.guess(10.*pars["composition"][elem]*pars["composition"][elem_scat], expr=f"n_{elem_scat:s}{elem:s}"))

#######################
### Execute Fitting ###
#######################

endloop = False
count = 0
while (not endloop) and (count<=2):
    starttime = datetime.datetime.now()
    print("start : ", starttime)

    out = larch.xafs.feffit(exafspars, dsets.flatten(), xtol=10**-xtol, ftol=10**-ftol, gtol=10**-gtol)

    report = larch.xafs.feffit_report(out)

    endtime = datetime.datetime.now()
    print("end   : ", endtime)
    print("time : ", endtime-starttime)

    count += 1
    if out.errorbars and count>1:
        endloop = True
    else:
        out_trial = f"{export_dir:s}/trial{count:d}"
        os.makedirs(f"{export_dir:s}/trial{count:d}", exist_ok=True)
        tools.export_params(exafspars, out, sname, export_path=out_trial)
        with open(f"{out_trial:s}/{sname:s}_fit.txt", mode="w") as f:
            f.write(f"start: {starttime}\nend: {endtime}\ndelta: {endtime-starttime}\n\n\n")
            f.write(report)
        f.close()
        tools.export_rss(out, sname, elems, temps=temps, export_path=out_trial)
        tools.export_n(out, sname, elems, pars, export_path=out_trial)
        tools.fit_eins(out_trial, sname, elems, elemcolors, np.array([tempsdict[t] for t in temps]))
        tools.fix_n(exafspars, out, dsets.flatten(), ftol)

######################
### Export Results ###
######################

for i in range(dsets.shape[0]):
    for j in range(dsets.shape[1]):
        dname = dnames[i,j]
        dset = dsets[i,j]
        elem = elems[i]
        outpath = f"{export_dir:s}/fit/{dname:s}_fit.png"
        tools.export_result(dset, pars.loc[elem], save_flag=True, export_path=outpath)

outpath = f"{export_dir:s}/{sname:s}_fit.txt"
with open(outpath, mode="w") as f:
    f.write(f"start: {starttime}\nend: {endtime}\ndelta: {endtime-starttime}\n\n\n")
    f.write(report)
f.close()

tools.export_rss(out, sname, elems, temps=temps, export_path=export_dir)
tools.export_n(out, sname, elems, pars, export_path=export_dir)
tools.export_params(exafspars, out, sname, export_path=export_dir)

tools.fit_eins(export_dir, sname, elems, elemcolors, np.array([tempsdict[t] for t in temps]))
