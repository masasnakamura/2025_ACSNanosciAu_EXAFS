import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
import larch
import xraydb
from scipy.constants import atomic_mass, hbar, Boltzmann
from lmfit import Parameters, minimize, fit_report

###matplotlib.use('QtAgg')

HBAR = hbar
AMU = atomic_mass
KB = Boltzmann

CENTI2INCHI = 1.0/2.54

plt.rcParams['font.size'] = 8
plt.rcParams['font.family']= "sans-serif"
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams['font.sans-serif'] = ["Arial"]
plt.rcParams["xtick.major.size"] = 6.0
plt.rcParams["ytick.major.size"] = 6.0
plt.rcParams["xtick.minor.size"] = 3.0
plt.rcParams["ytick.minor.size"] = 3.0

pd.set_option('display.max_rows', 900)
pd.set_option('display.max_columns', 900)

def fig_norm(dat: larch.Group, pars:dict, save_flag=False, export_path=None, show_flag=False):
    export_path = os.path.splitext(export_path)[0]
    
    export = pd.DataFrame({
        "energy"        : pd.Series(dat.energy),
        "mu"            : pd.Series(dat.mu),
        "pre_edge"      : pd.Series(dat.pre_edge),
        "post_edge"     : pd.Series(dat.post_edge),
        "bkg"           : pd.Series(dat.bkg),
        "norm"          : pd.Series(dat.norm),
        "flat"          : pd.Series(dat.flat),
        "k"             : pd.Series(dat.k),
        "chik"          : pd.Series(dat.chi),
        "k3chik"        : pd.Series(dat.k**3*dat.chi),
        "kwin"          : pd.Series(dat.kwin),
        "r"             : pd.Series(dat.r),
        "chir_mag"      : pd.Series(dat.chir_mag),
        "chir_re"       : pd.Series(dat.chir_re),
        "chir_im"       : pd.Series(dat.chir_im),
        "rwin"          : pd.Series(dat.rwin),
        "q"             : pd.Series(dat.q),
        "chiq_mag"      : pd.Series(dat.chiq_mag),
        "chiq_re"       : pd.Series(dat.chiq_re),
        "chiq_im"       : pd.Series(dat.chiq_im),
    })

    if export_path!=None:
        export_csvpath = Rf"{export_path:s}.csv"
        export.to_csv(export_csvpath, index=False)

    fig = plt.figure(figsize=(15*CENTI2INCHI,12*CENTI2INCHI))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    chik3_max = 0
    chir_mag_max = 0
    chiq_re_max = 0

    ax1.plot(dat.energy, dat.post_edge, linewidth=1, color="blueviolet")
    ax1.plot(dat.energy, dat.pre_edge, linewidth=1, color="olivedrab")
    ax1.plot(dat.energy, dat.bkg, linewidth=1, color="red")
    ax1.plot(dat.energy, dat.mu, linewidth=1, color="black")
    ax1.axvspan(dat.e0 + pars["pre_min"], dat.e0 + pars["pre_max"], color="coral", alpha=0.2, zorder=-100)
    ax1.axvspan(dat.e0 + pars["post_min"], dat.e0 + pars["post_max"], color="coral", alpha=0.2, zorder=-100)
    ax1.axvline(dat.e0, linewidth=1, color="lightseagreen", zorder=-100)

    ax2.plot(dat.k, dat.k**3*dat.chi, marker=".", markersize=3, linewidth=1, color="black")
    if chik3_max < np.max(np.absolute(dat.k**3*dat.chi)):
        chik3_max = np.max(np.absolute(dat.k**3*dat.chi))

    ax3.plot(dat.r, dat.chir_mag, marker=".", markersize=3, linewidth=1, color="black")
    if chir_mag_max < np.max(dat.chir_mag):
        chir_mag_max = np.max(dat.chir_mag)

    ax4.plot(dat.q, dat.chiq_re, marker=".", markersize=3, linewidth=1, color="black")
    if chiq_re_max < np.max(np.absolute(dat.chiq_re)):
        chiq_re_max = np.max(np.absolute(dat.chiq_re))

    ax2.plot(dat.k, dat.kwin*chik3_max*1.1, linewidth=1, color="navy")
    ax3.plot(dat.r, dat.rwin*chir_mag_max*1.1, linewidth=1, color="navy")
    ax4.plot(dat.k, dat.kwin*chiq_re_max*1.1, linewidth=1, color="navy")
    ax2.set_xlim(0,18)
    ax3.set_xlim(0,6)
    ax4.set_xlim(0,18)
    ax1.set_xlabel("$E$ / eV")
    ax2.set_xlabel("$k$ / $\mathrm{\AA}^{-1}$")
    ax3.set_xlabel("$R$ / $\mathrm{\AA}$")
    ax4.set_xlabel("$q$ / $\mathrm{\AA}^{-1}$")
    ax1.set_ylabel("$\mu t$")
    ax2.set_ylabel("$k^3 \chi (k)$ / $\mathrm{\AA}^{-3}$")
    ax3.set_ylabel("$|\chi (R)|$ / $\mathrm{\AA}^{-4}$")
    ax4.set_ylabel("$\mathrm{Re} [\chi (q)]$ / $\mathrm{\AA}^{-3}$")

    fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.15, wspace=0.5, hspace=0.3)

    if save_flag and export_path!=None:
        fig.savefig(Rf"{export_path:s}.png", transparent=True, dpi=600)
    if show_flag:
        plt.show()
    plt.clf()
    plt.close()

def fig_temps(exportdir: str, sname:str, elem: str, temps: list, tempcolors: dict, save_flag=False, show_flag=False):
    fig = plt.figure(figsize=(15*CENTI2INCHI,12*CENTI2INCHI))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    chik3_max = 0
    chir_mag_max = 0
    chiq_re_max = 0

    for temp in temps:
        dat = pd.read_csv(f"{exportdir:s}/dat/{sname:s}_{temp:s}_{elem:s}.csv")

        ax1.plot(dat["k"], dat["k"]**3*dat["chik"], marker=".", markersize=3, linewidth=1, color=tempcolors[temp])
        if chik3_max < np.max(np.absolute(dat["k"]**3*dat["chik"])):
            chik3_max = np.max(np.absolute(dat["k"]**3*dat["chik"]))

        ax2.plot(dat["r"], dat["chir_mag"], marker=".", markersize=3, linewidth=1, color=tempcolors[temp])
        if chir_mag_max < np.max(dat["chir_mag"]):
            chir_mag_max = np.max(dat["chir_mag"])

        ax3.plot(dat["q"], dat["chiq_re"], marker=".", markersize=3, linewidth=1, color=tempcolors[temp])
        if chiq_re_max < np.max(np.absolute(dat["chiq_re"])):
            chiq_re_max = np.max(np.absolute(dat["chiq_re"]))

    ax1.plot(dat["k"], dat["kwin"]*chik3_max*1.1, linewidth=1, color="navy")
    ax2.plot(dat["r"], dat["rwin"]*chir_mag_max*1.1, linewidth=1, color="navy")
    ax3.plot(dat["k"], dat["kwin"]*chiq_re_max*1.1, linewidth=1, color="navy")
    ax1.set_xlim(0,18)
    ax2.set_xlim(0,6)
    ax3.set_xlim(0,18)
    ax1.set_xlabel("$k$ / $\mathrm{\AA}^{-1}$")
    ax2.set_xlabel("$R$ / $\mathrm{\AA}$")
    ax3.set_xlabel("$q$ / $\mathrm{\AA}^{-1}$")
    ax1.set_ylabel("$k^3 \chi (k)$ / $\mathrm{\AA}^{-3}$")
    ax2.set_ylabel("$|\chi (R)|$ / $\mathrm{\AA}^{-4}$")
    ax3.set_ylabel("$\mathrm{Re} [\chi (q)]$ / $\mathrm{\AA}^{-3}$")

    fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)

    if save_flag:
        fig.savefig(f"{exportdir:s}/dat/{sname:s}_{elem:s}.png", transparent=True, dpi=600)
    if show_flag:
        plt.show()
    plt.clf()
    plt.close()

def export_result(dset:larch.Group, pars:dict, save_flag=False, export_path=None, show_flag=False):
    export_path = os.path.splitext(export_path)[0]

    figft = plt.figure(figsize=(15*CENTI2INCHI,12*CENTI2INCHI))
    larch.xafs.xftr(dset.data, rmin=pars["rmin_ftr"], rmax=pars["rmax_ftr"], dr=pars["dr_ftr"], window="hanning")
    larch.xafs.xftr(dset.model, rmin=pars["rmin_ftr"], rmax=pars["rmax_ftr"], dr=pars["dr_ftr"], window="hanning")
    export = pd.DataFrame({
        "mea_k"         : pd.Series(dset.data.k),
        "mea_chik"      : pd.Series(dset.data.chi),
        "mea_k3chik"    : pd.Series(dset.data.k**3*dset.data.chi),
        "mea_r"         : pd.Series(dset.data.r),
        "mea_chir_mag"  : pd.Series(dset.data.chir_mag),
        "mea_chir_re"   : pd.Series(dset.data.chir_re),
        "mea_chir_im"   : pd.Series(dset.data.chir_im),
        "mea_q"         : pd.Series(dset.data.q),
        "mea_chiq_mag"  : pd.Series(dset.data.chiq_mag),
        "mea_chiq_re"   : pd.Series(dset.data.chiq_re),
        "mea_chiq_im"   : pd.Series(dset.data.chiq_im),
        "fit_k"         : pd.Series(dset.model.k),
        "fit_chik"      : pd.Series(dset.model.chi),
        "fit_k3chik"    : pd.Series(dset.model.k**3*dset.model.chi),
        "fit_r"         : pd.Series(dset.model.r),
        "fit_chir_mag"  : pd.Series(dset.model.chir_mag),
        "fit_chir_re"   : pd.Series(dset.model.chir_re),
        "fit_chir_im"   : pd.Series(dset.model.chir_im),
        "fit_q"         : pd.Series(dset.model.q),
        "fit_chiq_mag"  : pd.Series(dset.model.chiq_mag),
        "fit_chiq_re"   : pd.Series(dset.model.chiq_re),
        "fit_chiq_im"   : pd.Series(dset.model.chiq_im),
        "win_k"         : pd.Series(dset.data.kwin),
        "win_r"         : pd.Series(dset.data.rwin)
    })

    if export_path!=None:
        export_csvpath = Rf"{export_path:s}.csv"
        export.to_csv(export_csvpath, index=False)

    axft2 = figft.add_subplot(2,2,1)
    axft3 = figft.add_subplot(2,2,3)
    axft4 = figft.add_subplot(2,2,2)
    axft5 = figft.add_subplot(2,2,4)

    axft2.plot(dset.data.r, dset.data.chir_mag, color="black", linewidth=0.5)
    axft2.plot(dset.data.r, 1.1*np.max(dset.data.chir_mag)*dset.data.rwin, color="gray")
    axft2.plot(dset.model.r, dset.model.chir_mag, color="red")
    axft2.scatter(dset.data.r, dset.data.chir_mag, s=5, color="black")
    axft2.set_xlim((0,6))

    axft3.plot(dset.data.r, dset.data.chir_re, color="black", linewidth=0.5)
    axft3.plot(dset.data.r, 1.1*np.max(dset.data.chir_re)*dset.data.rwin, color="gray")
    axft3.plot(dset.model.r, dset.model.chir_re, color="red")
    axft3.scatter(dset.data.r, dset.data.chir_re, s=5, color="black")
    axft3.set_xlim((0,6))

    axft4.plot(dset.data.k, dset.data.k**3*dset.data.chi, color="black", linewidth=0.5)
    axft4.plot(dset.data.k, 1.1*np.max(dset.data.k**3*dset.data.chi)*dset.data.kwin, color="gray")
    axft4.plot(dset.model.k, dset.model.k**3*dset.model.chi, color="red")
    axft4.scatter(dset.data.k, dset.data.k**3*dset.data.chi, s=5, color="black")
    axft4.set_xlim((0,np.max(dset.data.k)))

    axft5.plot(dset.data.q, dset.data.chiq_re, color="black", linewidth=0.5)
    axft5.plot(dset.data.k, 1.1*np.max(dset.model.chiq_re)*dset.data.kwin, color="gray")
    axft5.plot(dset.model.q, dset.model.chiq_re, color="red")
    axft5.scatter(dset.data.q, dset.data.chiq_re, s=5, color="black")
    axft5.set_xlim((0,np.max(dset.data.k)))

    axft3.set_xlabel("$R$ / $\mathrm{\AA}$")
    axft5.set_xlabel("$k,q$ / $\mathrm{\AA}^{-1}$")
    axft2.set_ylabel("$|\chi (R)|$ / $\mathrm{\AA}^{-4}$")
    axft3.set_ylabel("$\mathrm{Re} [\chi (R)]$ / $\mathrm{\AA}^{-4}$")
    axft4.set_ylabel("$k^3 \chi (k)$ / $\mathrm{\AA}^{-3}$")
    axft5.set_ylabel("$\mathrm{Re} [\chi (q)]$ / $\mathrm{\AA}^{-3}$")

    axft2.set_xlim(0,6)
    axft3.set_xlim(0,6)
    axft4.set_xlim(0,18)
    axft5.set_xlim(0,18)

    figft.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)

    if save_flag and export_path!=None:
        export_figpath = Rf"{export_path:s}.png"
        figft.savefig(export_figpath, transparent=True, dpi=600)
    if show_flag:
        plt.show()
    plt.clf()
    plt.close()

def export_rss(result, sname, elems, temps=None, export_path=None):
    if export_path!=None:
        export_path = os.path.splitext(export_path)[0]
        r_mean = np.empty((len(elems), len(temps)))
        r_err = np.empty((len(elems), len(temps)))
        ss_mean = np.empty((len(elems), len(temps)))
        ss_err = np.empty((len(elems), len(temps)))
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        df3 = pd.DataFrame()
        df4 = pd.DataFrame()
        for i,elem in enumerate(elems):
            for j,temp in enumerate(temps):
                r_param = getattr(result.paramgroup, f"r_{elem:s}_{temp:s}")
                r_mean[i,j] = r_param.value
                r_err[i,j] = r_param.stderr
                ss_param = getattr(result.paramgroup, f"ss_{elem:s}_{temp:s}")
                ss_mean[i,j] = ss_param.value
                ss_err[i,j] = ss_param.stderr
        for i,elem in enumerate(elems):
            df1i = pd.DataFrame(
                {
                    f"{elem:s}": r_mean[i,:],
                }
            )
            df2i = pd.DataFrame(
                {
                    f"{elem:s}": r_err[i,:],
                }
            )
            df3i = pd.DataFrame(
                {
                    f"{elem:s}": ss_mean[i,:],
                }
            )
            df4i = pd.DataFrame(
                {
                    f"{elem:s}": ss_err[i,:],
                }
            )
            df1 = pd.concat([df1, df1i], axis=1)
            df1 = df1.copy()
            df2 = pd.concat([df2, df2i], axis=1)
            df2 = df2.copy()
            df3 = pd.concat([df3, df3i], axis=1)
            df3 = df3.copy()
            df4 = pd.concat([df4, df4i], axis=1)
            df4 = df4.copy()
        df1.index = temps
        df2.index = temps
        df3.index = temps
        df4.index = temps

        df1.to_csv(f"{export_path:s}/{sname:s}_r_val.csv")
        df2.to_csv(f"{export_path:s}/{sname:s}_r_err.csv")
        df3.to_csv(f"{export_path:s}/{sname:s}_ss_val.csv")
        df4.to_csv(f"{export_path:s}/{sname:s}_ss_err.csv")

def export_n(result, sname, elems, pars, export_path=None):
    if export_path!=None:
        n_mean = np.empty((len(elems), len(elems)))
        n_err = np.empty((len(elems), len(elems)))
        for i,elem1 in enumerate(elems):
            for j,elem2 in enumerate(elems):
                n_param = getattr(result.paramgroup, f"n_{elem1:s}{elem2:s}")
                n_mean[i,j] = n_param.value/pars["composition"][elem1]
                try:
                    n_err[i,j] = n_param.stderr/pars["composition"][elem1]
                except TypeError:
                    n_err[i,j] = n_param.stderr
        df_mean = pd.DataFrame(n_mean, index=elems, columns=elems)
        df_err = pd.DataFrame(n_err, index=elems, columns=elems)
        df_mean["sum"] = df_mean.sum(axis=1)

        export_meanpath = f"{export_path:s}/{sname:s}_n_val.csv"
        export_errpath = f"{export_path:s}/{sname:s}_n_err.csv"
        df_mean.to_csv(export_meanpath)
        df_err.to_csv(export_errpath)

def export_log(python_path, par_norm, tol=None, fit_init=None, log_path=None):
    if log_path!=None:
        with open(python_path, encoding="UTF-8") as f:
            python_code = f.read()
        f.close()

        with open(log_path, mode="w") as f:
            f.write("### PARAMETERS ###\n")
            f.write(str(par_norm))
            f.write("\n\n\n")

            f.write("### FITTING ###\n")
            if fit_init!=None:
                f.write(f"path = {str(fit_init):s}\n")
            if tol!=None:
                f.write(f"tol = {tol:d}\n")

            f.write("### CODE USED ###\n")
            f.write(python_code)
            f.write("\n\n\n")
        f.close()

def norm_frac_edge(dat, pars, e0=None, elem=None, edge=None, frac=0.5):
    count=0
    cri=10
    eold=0
    if elem!=None and edge!=None:
        e0 = xraydb.xray_edge(elem, edge, energy_only=True)
    larch.xafs.pre_edge(dat.energy, dat.mu, e0=e0, group=dat, pre1=pars["pre_min"], pre2=pars["pre_max"], norm1=pars["post_min"], norm2=pars["post_max"], nnorm=pars["nnorm"])
    while cri > 0.001 and count < 100:
        for k, e in enumerate(dat.energy):
            if dat.flat[k] < frac < dat.flat[k+1]:
                enew = dat.energy[k] + (dat.energy[k+1] - dat.energy[k]) * (frac-dat.flat[k]) / (dat.flat[k+1]-dat.flat[k])
                break
        cri = np.abs(enew - eold)
        eold = enew
        larch.xafs.pre_edge(dat.energy, dat.mu, e0=enew, group=dat, pre1=pars["pre_min"], pre2=pars["pre_max"], norm1=pars["post_min"], norm2=pars["post_max"], nnorm=pars["nnorm"])
        count+=1

    dat.e0 = enew
    return

def fit_eins(outdir, sname, elems, elemcolors, temps: np.ndarray):
    def eins(pars, x):
        return ((10**20)*(HBAR**2.)/(2.*pars["m"]*AMU*KB*pars["thetae"])/np.tanh(pars["thetae"]/(2.*x))) + pars["ss_stat"]

    def objective(pars, x, y, y_err):
        return (y - eins(pars, x))/y_err

    outpath = f"{outdir:s}/{sname:s}_ss_fit.txt"

    with open(outpath, mode="w") as f:
        f.write("")
    f.close()

    fig, axes = plt.subplots(2, 3, figsize=(15*CENTI2INCHI,12*CENTI2INCHI))

    datval = pd.read_csv(f"{outdir:s}/{sname:s}_ss_val.csv", header=0, index_col=0)
    datstd = pd.read_csv(f"{outdir:s}/{sname:s}_ss_err.csv", header=0, index_col=0)

    df = pd.DataFrame(columns=["val", "std"])

    for i,elem in enumerate(elems):
        params = Parameters()
        params.add("thetae", value=200., vary=True)
        params.add("ss_stat", value=0., vary=True)
        params.add("m", value=xraydb.atomic_mass(elem), vary=False)

        stdflag = True

        try:
            out = minimize(objective, params, args=(temps, datval[elem], datstd[elem]), scale_covar=False)
        except ValueError:
            out = minimize(objective, params, args=(temps, datval[elem], 1.), scale_covar=False)
            stdflag = False

        report = fit_report(out)

        with open(outpath, mode="a") as f:
            f.write(f"### {elem:s} ###\n")
            f.write(report)
            f.write("\n\n\n")
        f.close()

        if stdflag:
            axes[i//3, i%3].errorbar(temps, datval[elem], datstd[elem], fmt="o", color=elemcolors[elem])
        else:
            axes[i//3, i%3].scatter(temps, datval[elem], color=elemcolors[elem])
        axes[i//3, i%3].plot(np.arange(1.,401.,1.), eins(out.params, np.arange(1., 401., 1.)), color=elemcolors[elem], alpha=0.5, zorder=-100)
        axes[i//3, i%3].text(0.05, 0.95, r"$\theta_\mathrm{E}$ = %.0f(%.0f) K" % (out.params["thetae"].value, out.params["thetae"].stderr), ha="left", va="top", transform=axes[i//3, i%3].transAxes)
        axes[i//3, i%3].set_xlim(-5,(np.max(temps))//100*100+105)
        if i%3==0:
            axes[i//3, i%3].set_ylabel("$\sigma^2$ / $\mathrm{\AA}^2$")
        if i//3==1:
            axes[i//3, i%3].set_xlabel("$T$ / K")

        df.loc[elem] = [out.params["thetae"].value, out.params["thetae"].stderr]
    
    if len(elems)==5:
        axes[1,2].axis("off")

    fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)
    fig.savefig(f"{outdir:s}/{sname:s}_ss_fit.png", dpi=600, transparent=True)
    plt.clf()
    plt.close()

    df.to_csv(f"{outdir:s}/{sname:s}_ss_fit.csv")


def fix_n(pars, result, datasets, ftol):
    def resid(params, datasets):
        return np.concatenate([d._residual(params) for d in datasets])
    for p in pars.__dir__():
        getattr(pars,p).value = result.params[p].value
        if (p[:2]=="n_"):
            pars2 = copy.deepcopy(pars)
            getattr(pars2,p).value = 0.
            getattr(pars2,p).vary = False
            cost1 = np.sum(resid(pars, datasets)**2.0)
            cost2 = np.sum(resid(pars2, datasets)**2.0)
            if np.abs(cost2-cost1)/cost1 < 10**-(ftol+1):
                getattr(pars,p).value = 0.
                getattr(pars,p).vary = False

def export_params(pars, result, sname, export_path):
    df = pd.DataFrame(columns=["value", "stderr"])
    for p in pars.__dir__():
        df.loc[p] = [result.params[p].value, result.params[p].stderr]
    df.to_csv(f"{export_path:s}/{sname:s}_params.csv")