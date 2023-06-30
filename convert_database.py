from sage.all import *
import sys
import os
from pyhdme import modular_igusa_from_igusa_clebsch

from collections import defaultdict
# st_to_geomendR = {'USp(4)': 'R', 'SU(2)xSU(2)': 'RxR', 'N(SU(2)xSU(2))': 'RxR', 'U(1)xSU(2)': 'CxR', 'N(U(1)xSU(2))': 'CxR', 'F': 'CxC', 'F_{ab}': 'CxC', 'F_a': 'CxC', 'F_{ac}': 'CxC', 'F_{a,b}': 'CxC', 'E_1': 'M_2(R)', 'E_2': 'M_2(R)', 'J(E_1)': 'M_2(R)', 'E_3': 'M_2(R)', 'E_4': 'M_2(R)', 'J(E_2)': 'M_2(R)', 'J(E_3)': 'M_2(R)', 'E_6': 'M_2(R)', 'J(E_4)': 'M_2(R)', 'J(E_6)': 'M_2(R)', 'C_1': 'M_2(C)', 'J(C_1)': 'M_2(C)', 'C_2': 'M_2(C)', 'C_{2,1}': 'M_2(C)', 'C_3': 'M_2(C)', 'C_{4,1}': 'M_2(C)', 'C_4': 'M_2(C)', 'D_2': 'M_2(C)', 'J(C_2)': 'M_2(C)', 'D_{2,1}': 'M_2(C)', 'D_3': 'M_2(C)', 'D_{3,2}': 'M_2(C)', 'J(C_3)': 'M_2(C)', 'C_{6,1}': 'M_2(C)', 'C_6': 'M_2(C)', 'J(C_4)': 'M_2(C)', 'D_{4,1}': 'M_2(C)', 'D_4': 'M_2(C)', 'D_{4,2}': 'M_2(C)', 'J(D_2)': 'M_2(C)', 'T': 'M_2(C)', 'J(D_3)': 'M_2(C)', 'D_{6,1}': 'M_2(C)', 'D_6': 'M_2(C)', 'D_{6,2}': 'M_2(C)', 'J(C_6)': 'M_2(C)', 'J(D_4)': 'M_2(C)', 'O_1': 'M_2(C)', 'O': 'M_2(C)', 'J(T)': 'M_2(C)', 'J(D_6)': 'M_2(C)', 'J(O)': 'M_2(C)'}


def coefficients(eqn_str):
    R = PolynomialRing(Integers(), 'x')
    if 'x' in eqn_str:
        coeffs = list(map(list, map(R, eqn_str.strip('[]').split(","))))
    else:
        coeffs = [[] if not elt else list(map(int, elt.split(','))) for elt in  eqn_str.lstrip('[').rstrip(']').split("],[")]
    return coeffs

def modular_inv(coeffs):
    from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants
    R = PolynomialRing(Rationals(), 'x')
    f, g = map(R, coeffs)
    return modular_igusa_from_igusa_clebsch(igusa_clebsch_invariants(4*f + g**2))


from multiprocessing import cpu_count
@parallel(ncpus=cpu_count())
def doline(lines):
    res = []
    for line in lines:
        _, cond, Lhash, eqn_str, _, _, ST =  line.strip().split(':')
        if ST != 'USp(4)':
            continue
        cond, Lhash = map(int, (cond, Lhash))
        coeffs = coefficients(eqn_str)
        inv = modular_inv(coeffs)
        res.append((cond, Lhash, inv, coeffs))
    return res


def load_drew(filename):
# edit accordingly
    sys.path.append(os.path.join(os.environ['HOME'],"projects/LMFDB/lmfdb"))
    from lmfdb import db
    res = defaultdict(list)
    names = {int(elt['Lhash']) : elt['class'] for elt in db.g2c_curves.search({}, ['Lhash', 'class'])}

    print(f"Processing {filename}...")
    with open(filename) as F:
        lines = F.readlines()
        blocksoflines = [lines[i:len(lines):cpu_count()] for i in range(cpu_count())]
        assert sum(len(elt) for elt in blocksoflines)
        for _, rlines in doline(blocksoflines):
            for (cond, Lhash, inv, coeffs) in rlines:
                res[(cond, Lhash)].append([list(inv), coeffs])
    print(f"\rDone")
    for elt in res:
        res[elt].sort()
        res[elt] = [[ elt[0] for elt in res[elt]], [elt[1] for elt in res[elt]]]
    return dict(res), dict(names)


def generate_input(filename, database, names):
    print(f"Generating {filename}...")
    with open(filename, "w") as W:
        for k in sorted(database):
            v = database[k]
            cond, Lhash = k
            name = names.get(Lhash)
            row = [cond, name, Lhash] + v
            W.write(":".join(map(str, row)).replace(" ", "") + "\n")
    print(f"\rDone")


def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')

if __name__ == "__main__" and not is_interactive():
    if len(sys.argv) != 3:
        print(r"""
    USAGE:
         sage -python {cmd} <input filename> <output filename>

         where the input filename should be a downloaded version of https://math.mit.edu/~drew/st2e20.txt

    EXAMPLE
        """.format(cmd=sys.argv[0]) )
        sys.exit()

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    drewst4, names = load_drew(input_filename)
    generate_input(output_filename, drewst4, names)

