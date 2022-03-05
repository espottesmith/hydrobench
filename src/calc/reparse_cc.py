# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks.core.launchpad import LaunchPad
from fireworks import Workflow
from pymatgen.core.structure import Molecule
import copy
import os

from pymongo import UpdateOne
from pymongo.errors import OperationFailure

from monty.json import MSONable, jsanitize
from monty.serialization import loadfn, dumpfn

from fireworks import Firework, Workflow, LaunchPad

from atomate.common.powerups import add_tags
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.fireworks.core import SinglePointFW
from atomate.qchem.database import QChemCalcDb
from atomate.qchem.drones import QChemDrone

db = QChemCalcDb.from_db_file("/Users/ewcss/config/fireworks/ewcss_db.json")
lp = LaunchPad.from_file("/Users/ewcss/config/fireworks/ewcss_launchpad.yaml")

docs = list()
for calc in list(db.db["tasks"].find({"tags.set": "20211016_sp_cc",
                                      "dir_name": {"$regex": "eaglefs"},
                                      },
                                     {"task_id": 1, "task_label": 1, "dir_name": 1, "tags": 1})):

    t = [QChemToDb(calc_dir=calc["dir_name"],
                   db_file="/home/ewcss/config/atomate/ewcss_db.json",
                   additional_fields={"task_label": calc["task_label"]})]
    fw = Firework(t, name=calc["task_label"])
    wf = Workflow([fw], name=calc["task_label"])
    wf = add_tags(wf, {"class": calc["tags"]["class"], "set": "20211016_sp_cc"})
    print(wf)
    lp.add_wf(wf)