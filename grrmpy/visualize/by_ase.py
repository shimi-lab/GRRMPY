import ase
import ase.visualize
from ase import Atoms
import ipywidgets as widgets
from ipywidgets import interactive
from IPython.display import display

def view(atoms):
    """再生ボタンのない構造を表示"""
    if type(atoms)==list:
        _view_atomslist(atoms)
    elif type(atoms)==Atoms:
        """Atomsで与えられた場合"""
        _view_atoms(atoms)
    else:
        raise TypeError("引数が間違っています")
    
def _view_atoms(atoms:Atoms):
        v = ase.visualize.view(atoms, viewer='ngl')
        v.view.add_representation("ball+stick")
        display(v)  
        
def _view_atomslist(atomslist:list):
        """Atomsリストで与えられた場合""" 
        def disp_structure(i):
            if atomslist[i]:
                v = ase.visualize.view(atomslist[i], viewer='ngl')
            else:
                v = ase.visualize.view(Atoms(), viewer='ngl')
            v.view.add_representation("ball+stick")
            display(v)
                       
        item = widgets.BoundedIntText(
            value=0,
            min=0,
            max=len(atomslist)-1,
            step=1,
            description=' ',
            disabled=False
        )
        slider = widgets.IntSlider(
            value=0,
            min=0,
            max=len(atomslist)-1,
            step=1,
            description=' ',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='d'
        )        
        mylink = widgets.link((item, 'value'), (slider, 'value'))
        w = interactive(disp_structure, i=slider)
        display(item)
        display(w) 
    
def view_images(imagess):
    """再生ボタン付きの構造を表示"""
    if type(imagess)==list:
        if type(imagess[0])== Atoms:
            """再生ボタン付き"""
            _viewimage_1datomslist(imagess)
        else:
            _viewimage_2datomslist(imagess)
    else:
        raise TypeError("引数が間違っています")

def _viewimage_1datomslist(images):
    v = ase.visualize.view(images, viewer='ngl')
    v.view.add_representation("ball+stick")
    display(v)

def _viewimage_2datomslist(imagess):
        def disp_structure(i):
            if imagess[i]:
                v = ase.visualize.view(imagess[i], viewer='ngl')
                v.view.add_representation("ball+stick")
            else:
                v = "データなし"
            display(v)
        item = widgets.BoundedIntText(
            value=0,
            min=0,
            max=len(imagess)-1,
            step=1,
            description=' ',
            disabled=False
        )
        slider = widgets.IntSlider(
            value=0,
            min=0,
            max=len(imagess)-1,
            step=1,
            description=' ',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='d'
        )        
        mylink = widgets.link((item, 'value'), (slider, 'value'))
        w = interactive(disp_structure, i=slider)
        display(item)
        display(w)