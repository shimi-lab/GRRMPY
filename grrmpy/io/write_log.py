from ase.units import Hartree
# USER
from grrmpy.conv.atoms2coordination import atoms2coordination
from grrmpy.calculator import pfp_calculator

def write_eqlog(file,atoms_list,indices=None,calc_func=pfp_calculator):
    with open(file,"w") as f:
        f.write("List of Equilibrium Structures\n\n")
        for i,atoms in enumerate(atoms_list):
            f.write(f"# Geometry of EQ {i}, SYMMETRY = -- \n")
            f.write(atoms2coordination(atoms,indices))
            energy = _get_energy(atoms,calc_func)
            f.write(f"Energy    =  {energy} ( {energy} :    0.000000000000)\n\n")
            
def write_tslog(file,atoms_list,connections_list,ts=True,indices=None,calc_func=pfp_calculator):
    first_txt = "List of Transition Structures\n\n" if ts else "List of Path Top (Approximate TS) Structures\n\n"
    with open(file,"w") as f:
        f.write(first_txt)
        for i,(atoms,(ini,fin)) in enumerate(zip(atoms_list,connections_list)):
            f.write(f"# Geometry of TS {i}, SYMMETRY = --  \n")
            f.write(atoms2coordination(atoms,indices))
            energy = _get_energy(atoms,calc_func)
            f.write(f"Energy    =  {energy} ( {energy} :    0.000000000000)\n")
            f.write(f"CONNECTION : {ini} - {fin}\n\n")
    
def _get_energy(atoms,calc_func=pfp_calculator):
    try:
        return atoms.get_potential_energy()/Hartree
    except:
        atoms.calc = calc_func()
        return atoms.get_potential_energy()/Hartree