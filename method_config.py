"""
Gaussian method for the ArX calculations. Read by 03_make_gaussian_inputs.py.

Route lines and title copied from g16_inpgen_arbrlib.py (default, arbr_lib = True), the script that wrote the
ArBr library inputs, and checked against an actual ArBr input (arbr15163_conf-1.com). Link0 (16 cores,
16 GB) is taken from that input; the script's own default was 24 cores / 64 GB:

  link 1: #N empiricaldispersion=GD3BJ OPT FREQ=NORAMAN Integral(Ultrafine) B3LYP 6-31G(d,p)
  link 2: #N int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check M062X def2TZVP
  link 3: #N int=(grid=ultrafine) nmr=giao guess=read geom=check M062X def2TZVP

"{basis_opt}" / "{basis_sp}" stand where the basis keyword goes. The script fills them in:
  ArCl -> 6-31G(d,p) and def2TZVP: identical to the ArBr inputs (both basis sets cover Cl).
  ArI  -> GenECP with the basis sets in HEAVY below, because 6-31G(d,p) is defined only for H-Kr.
"""

DRAFT = False

LINK0 = ["%NProcShared=16", "%Mem=16GB"]      # as in an actual ArBr input (arbr15163_conf-1.com); %chk is written after these

LINKS = [
    dict(route="#N empiricaldispersion=GD3BJ OPT FREQ=NORAMAN Integral(Ultrafine) B3LYP {basis_opt}",
         read_geometry=False, nbo=""),
    dict(route="#N int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check M062X {basis_sp}",
         read_geometry=True, nbo=""),
    dict(route="#N int=(grid=ultrafine) nmr=giao guess=read geom=check M062X {basis_sp}",
         read_geometry=True, nbo=""),
]

# opt + freq count as two jobs, so 2 + 1 + 1 = 4 (same as the old log checker's default)
EXPECTED_NORMAL_TERMINATIONS = 4

BASIS_OPT = "6-31G(d,p)"
BASIS_SP = "def2TZVP"

# Elements heavier than Kr (here: iodine). The ArBr library had none, so there is no ArBr precedent.
# !! Methods decision to confirm before running ArI !!
#   opt: 6-31G(d,p) on light atoms + def2-SVP with the def2 ECP on I
#   sp:  def2-TZVP on all atoms, with the def2 ECP on I requested explicitly
HEAVY = dict(opt_basis="def2SVP", sp_basis="def2TZVP", ecp="def2")
