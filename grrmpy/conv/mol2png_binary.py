from rdkit.Chem.Draw import rdMolDraw2D
import base64
from rdkit import rdBase

def mol2png_binary(mol, size=(250,250), kekulize=True, wedgeBonds=True, options=None, **kwargs):
    """molオブジェクトをpngのbynaryデータに変換しutf-8でエンコードした状態で返す
    
    Parameters:
    
    mol : Mol
        Molオブジェクト
    size : tuple
        画像のサイズ
        
    Retrun:
        str: バイナリデータをutf-8に変換したもの
    """
    data = _moltoimg(mol, size, kwargs.get('highlightAtoms', []), kwargs.get('legend', ''),
                highlightBonds=kwargs.get('highlightBonds', []), drawOptions=options,
                kekulize=kekulize, wedgeBonds=wedgeBonds, returnPNG=True)
    return base64.b64encode(data).decode('utf-8')

def _okToKekulizeMol(mol, kekulize):
    if kekulize:
        for bond in mol.GetBonds():
            if bond.GetIsAromatic() and bond.HasQuery():
                return False
        return True
    return kekulize

def _moltoimg(mol, sz, highlights, legend, returnPNG=False, drawOptions=None, **kwargs):
    try:
        with rdBase.BlockLogs():
            mol.GetAtomWithIdx(0).GetExplicitValence()
    except RuntimeError:
        mol.UpdatePropertyCache(False)

    kekulize = _okToKekulizeMol(mol, kwargs.get('kekulize', True))
    wedge = kwargs.get('wedgeBonds', True)

    try:
        with rdBase.BlockLogs():
            mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize, wedgeBonds=wedge)
    except ValueError:  # <- can happen on a kekulization failure
        mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False, wedgeBonds=wedge)

    d2d = rdMolDraw2D.MolDraw2DCairo(sz[0], sz[1])
    if drawOptions is not None:
        d2d.SetDrawOptions(drawOptions)
    if 'highlightColor' in kwargs and kwargs['highlightColor']:
        d2d.drawOptions().setHighlightColour(kwargs['highlightColor'])
    # we already prepared the molecule:
    d2d.drawOptions().prepareMolsBeforeDrawing = False
    bondHighlights = kwargs.get('highlightBonds', None)
    if bondHighlights is not None:
        d2d.DrawMolecule(mc, legend=legend or "", highlightAtoms=highlights or [],
                    highlightBonds=bondHighlights)
    else:
        d2d.DrawMolecule(mc, legend=legend or "", highlightAtoms=highlights or [])
    d2d.FinishDrawing()
    img = d2d.GetDrawingText()
    return img