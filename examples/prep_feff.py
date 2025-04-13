import os
import numpy as np
import xraydb
import larch
from larch import Interpreter
session = Interpreter()

eldict = dict(
    {
        "Ru":
        {
            "edge": "K",
            "r": 1.339
        },
        "Rh":
        {
            "edge": "K",
            "r": 1.345
        },
        "Pd":
        {
            "edge": "K",
            "r": 1.375
        },
        "Ir":
        {
            "edge": "L3",
            "r": 1.358
        },
        "Pt":
        {
            "edge": "L3",
            "r": 1.386
        },
        "Ga":
        {
            "edge": "K",
            "r": 1.354
        },
        "In":
        {
            "edge": "K",
            "r": 1.520
        },
        "Sn":
        {
            "edge": "K",
            "r": 1.476
        }
    }
)

templatep = f"./feffinp/template_fcc.inp"
with open(templatep, mode="r") as templatef:
    templates = templatef.read()
templatef.close()

elabss = eldict.keys()
elscats = eldict.keys()

hnumdict = {"K":"1", "L1":"2", "L2":"3", "L3":"4"}

for elabs in elabss:
    for elscat in elscats:
        zabs = xraydb.atomic_number(elabs)
        zscat = xraydb.atomic_number(elscat)
        edge = eldict[elabs]["edge"]
        hnum = hnumdict[edge]
        d = eldict[elabs]["r"]+eldict[elscat]["r"]

        exps = templates.replace("Zabs", str(zabs)).replace("Zsct", str(zscat)).replace("Elabs", elabs).replace("Elsct", elscat).replace("hnum", hnum).replace("R100", str(d/np.sqrt(2.))).replace("R110", str(d)).replace("R200", str(d*np.sqrt(2.))).replace("R211", str(d*np.sqrt(3.))).replace("R220", str(d*2.))

        expdir = f"./feffinp/{elabs:s}{edge:s}-{elscat:s}"
        os.makedirs(expdir, exist_ok=True)
        expfile = f"{expdir:s}/{elabs:s}{edge:s}-{elscat:s}.inp"
        with open(expfile, mode="w") as expf:
            expf.write(exps)

        larch.xafs.feff6l(folder=expdir, feffinp=f"{elabs:s}{edge:s}-{elscat:s}.inp")
