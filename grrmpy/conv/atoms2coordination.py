def atoms2coordination(atoms,indices=None, fmt='%22.15f'):
    if indices is None:
        indices = [i for i in range(len(atoms))]
    atoms = atoms[indices]
    text = ""
    for s, (x, y, z) in zip(atoms.symbols, atoms.positions):
        text += '%-2s %s %s %s\n' % (s, fmt % x, fmt % y, fmt % z)
    return text