from io import StringIO
import nglview as nv
from nglview import NGLWidget
from ase import Atoms
import numpy as np
from ase.io.proteindatabank import write_proteindatabank
from typing import Any,Dict,List,Optional,Union
from nglview.component import ComponentViewer
from math import sin,cos,radians


def generate_js_code(color_dict):
    """構造をViewerで可視化する際の原子の色を規定するJSコードを作成する.
    
    Parameters:
    
    color_dict: dict
        | 元素名をkey,hexカラーコードをvaluseとする辞書
        | ex) {'H':'0xffcccc','C':'0x804929'}
        | 先頭に0xを付けなければいけない事に注意

    Returns:
        str: JavaScriptコード
    """
    #Jsコード作成
    jscode = "this.atomColor = function (atom) {"
    for i,(element, color) in enumerate(color_dict.items()):
        if i == 0:
           jscode += 'if (atom.element == ' + f'"{element.upper()}"' + ') { return ' + f'{color}'
        else:
            jscode += '} else if (atom.element == ' + f'"{element.upper()}"'+') {return '+f'{color}'
    jscode += '}}'
    
    return jscode

    # grrmpy.io.read_elementiniを参照すると良い
    
def update_tooltip_atoms(view: NGLWidget, atoms: Atoms):
    """原子上にマウスを置いたときに,原子のindexと位置を表示する.
    """
    if atoms.get_pbc().any():
        _, Q = atoms.cell.standard_form()
    else:
        Q = np.eye(3)
    Q_str = str(Q.tolist())
    var_str = f"this._Q = {Q_str}"
    script_str = """
    var tooltip = document.createElement('div');
    Object.assign(tooltip.style, {
      display: 'none',
      position: 'fixed',
      zIndex: 10,
      pointerEvents: 'none',
      backgroundColor: 'rgba( 0, 0, 0, 0.6 )',
      color: 'lightgrey',
      padding: '8px',
      fontFamily: 'sans-serif'
    });
    document.body.appendChild(tooltip);

    var that = this;
    this.stage.mouseControls.remove('hoverPick');
    this.stage.signals.hovered.add(function (pickingProxy) {
      if (pickingProxy && (pickingProxy.atom || pickingProxy.bond)) {
        var atom = pickingProxy.atom || pickingProxy.closestBondAtom
        var mp = pickingProxy.mouse.position
        //tooltip.innerText = atom.element + ' i=' + atom.index + ' (' + atom.x.toFixed(2) +  ', ' + atom.y.toFixed(2) +  ', ' + atom.z.toFixed(2) + ')'
        //var pos = that._atoms_pos[atom.index]
        var Q = that._Q
        var pos_x = Q[0][0] * atom.x + Q[0][1] * atom.y + Q[0][2] * atom.z 
        var pos_y = Q[1][0] * atom.x + Q[1][1] * atom.y + Q[1][2] * atom.z
        var pos_z = Q[2][0] * atom.x + Q[2][1] * atom.y + Q[2][2] * atom.z
        tooltip.innerText = 'i=' + atom.index + ' ' + atom.element + ' (' + pos_x.toFixed(2) +  ', ' + pos_y.toFixed(2) +  ', ' + pos_z.toFixed(2) + ')'
        tooltip.style.bottom = window.innerHeight - mp.y + 3 + 'px'
        tooltip.style.left = mp.x + 3 + 'px'
        tooltip.style.display = 'block'
      } else {
        tooltip.style.display = 'none'
      }
    });
    this.stage.tooltip = tooltip;
    """
    # https://github.com/nglviewer/ngl/blob/bd4a31c72e007d170b6bae298a5f7c976070e173/src/stage/mouse-behavior.ts#L31-L33
    view._execute_js_code(var_str + script_str)
    
def get_struct(atoms: Atoms, ext="pdb", replace_resseq: bool = False) -> List[Dict]:
    """Convert from ase `atoms` to `struct` object for nglviewer"""
    if ext == "pdb":
        # Faster, by using `StringIO`.
        sio = StringIO("")
        write_proteindatabank(sio, atoms)
        struct_str = sio.getvalue()
    else:
        struct_str = nv.ASEStructure(atoms, ext=ext).get_structure_string()
    if replace_resseq:
        struct_str = _replace_resseq(struct_str)
    struct = [dict(data=struct_str, ext=ext)]
    return struct
  
def _replace_resseq(struct: str) -> str:
    """Overwrite residue sequence number as atom index

    Please refer http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    for the format definition.

    This method was used to show atom index in default tooltip by overwriting
    residue name as tentative workaround.
    However, residue group is used for calculating bond information,
    and this method will affect to bond calculation in ngl.

    Args:
        struct (str): pdb string before replace

    Returns:
        struct (str): pdb string after residue sequence number replaced by atom index.
    """
    lines = struct.split("\n")
    atom_index = 0
    for i, line in enumerate(lines):
        if line.startswith("ATOM  "):
            lines[i] = line[:22] + f"{atom_index % 10000:4}" + line[26:]
            atom_index += 1
    return "\n".join(lines)
  
def add_force_shape(
    atoms: Atoms,
    view: NGLWidget,
    force_scale: float = 0.5,
    force_color: Optional[List[int]] = None,
) -> ComponentViewer:
    if force_color is None:
        force_color = [1, 0, 0]  # Defaults to red color.
    # Add force components
    forces = atoms.get_forces()
    pos = atoms.positions
    if atoms.get_pbc().any():
        rcell, rot_t = atoms.cell.standard_form()
        rot = rot_t.T
        pos_frac = pos.dot(rot)
        force_frac = forces.dot(rot)
    else:
        pos_frac = pos
        force_frac = forces

    shapes = []
    for i in range(atoms.get_global_number_of_atoms()):
        pos1 = pos_frac[i]
        pos2 = pos1 + force_frac[i] * force_scale
        pos1_list = pos1.tolist()
        pos2_list = pos2.tolist()
        shapes.append(("arrow", pos1_list, pos2_list, force_color, 0.2))
    c = view._add_shape(shapes, name="Force")
    return c

def rotate_view(view, x=0, y=0, z=0, degrees=True):
    """Rotate view over the x, y and z angles.
    
    Rotate the given `view` object for the given angles `x`, 'y' and `z` over
    their respective axes.
    
    Parameters
    ----------
    view : :obj:`nglview.NGLWidget`
        View widget to operate on.
    x : float, optional
        Angle of rotation over x-axis, *i.e.*, roll (the default is 0, no rotation).
    y : float, optional
        Angle of rotation over y-axis, *i.e.*, pitch (the default is 0, no rotation).
    z : float, optional
        Angle of rotation over z-axis, *i.e.*, yaw (the default is 0, no rotation).
    degrees : bool, optional
        True if angle is in degrees, False if angle is in radians (the default is True).
    
    Notes
    -----
    See [wikipedia](https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations)
    for implementation details.
    
    Example
    -------
    >>> import nglview as nv
    >>> view = nv.show_structure_file(nv.datafiles.PDB)
    >>> rotate_view(view, 90, 45, 15)
    >>> view
    """
    # Compute rotation matrix
    R = rotmat(x, y, z, degrees)
    # Compute quaternion
    Q = rotm2quat(R)
    # Perform rotation
    view.control.rotate(Q)
    
def rotmat(x, y, z, degrees=True):
    """Create a general rotation matrix for a given set of angles.
    
    Create a general rotation matrix for a given set of angles `x`, `y` and `z`.
    
    Parameters
    ----------
    x : float
        Angle of rotation over x-axis (roll).
    y : float
        Angle of rotation over y-axis (pitch).
    z : float
        Angle of rotation over z-axis (yaw).
    degrees : bool, optional
        True if angle is in degrees, False if angle is in radians, default: True
        
    Returns
    -------
    R: 3x3 rotation matrix
    
    Notes
    -----
    See (wikipedia)[https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations]
    for details.
    
    Example
    -------
    >>> rotmat(90, 45, 30)
    array([[ 0.61237244,  0.61237244,  0.5       ],
           [ 0.35355339,  0.35355339, -0.8660254 ],
           [-0.70710678,  0.70710678,  0.        ]])
    """
    import numpy as np
    from numpy import sin, cos
    
    c, b, a = [float(x), float(y), float(z)]
    degrees = bool(degrees)
    
    # Convert to degrees if needed
    if degrees:
        a = a * np.pi / 180.0
        b = b * np.pi / 180.0
        c = c * np.pi / 180.0
    
    # Compute rotation matrix
    R = [ 
        [ 
            np.cos(a) * np.cos(b),
            np.cos(a) * np.sin(b) * np.sin(c) - np.sin(a) * np.cos(c),
            np.cos(a) * np.sin(b) * np.cos(c) + np.sin(a) * np.sin(c)
        ],
        [
            np.sin(a) * np.cos(b),
            np.sin(a) * np.sin(b) * np.sin(c) + np.cos(a) * np.cos(c),
            np.sin(a) * np.sin(b) * np.cos(c) - np.cos(a) * np.sin(c)
        ],
        [
            -np.sin(b),
            np.cos(b) * np.sin(c),
            np.cos(b) * np.cos(c)
        ]

    ]
    
    # Round the results to solve machine precision problems
    # This way 1 is 1 and 0 is 0, etc.
    R = np.around(np.array(R), decimals=15)
    
    return R

def rotm2quat(R):
    """Convert a 3x3 rotation matrix to a quaternion.
    
    Convert a 3x3 rotation matrix to a 4x1 quaternion, which in turn can be used
    to rotate an object.
    
    Parameters
    ----------
    R : array_like 
        Rotation matrix (3x3 matrix) to convert to a quaternion.
        
    Returns
    -------
    list 
        Quaternion (4x1 list) for rotation matrix `R`.
    
    Raises
    ------
    ValueError
        If `R` is not a 3x3 matrix.
    
    Notes
    -----
    See [this website](See: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/)
    for implementation details.
    
    Examples
    --------
    >>> R = [[ 1.,  0.,  0.], [ 0.,  0., -1.], [-0.,  1.,  0.]]
    >>> rotm2quat(R)
    [0.7071067811865476, 0.7071067811865475, 0.0, 0.0]
    
    """
    import numpy as np
    
    # Ensure R is a numpy array
    R = np.array(R)
    
    # Check that R is 3x3 matrix
    if R.shape != (3,3):
        raise ValueError('Invalid rotation matrix shape {}x{}'.format(*R.shape))
        
    tr = R[0,0] + R[1,1] + R[2,2]
    
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2 # S=4*qw
        qw = 0.25 * S
        qx = (R[2,1] - R[1,2]) / S
        qy = (R[0,2] - R[2,0]) / S
        qz = (R[1,0] - R[0,1]) / S
    elif (R[0,0] > R[1,1]) and (R[0,0] > R[2,2]):
        S = np.sqrt(1 + R[0,0] - R[1,1] - R[2,2]) * 2 # S=4*qx
        qw = (R[2,1] - R[1,2]) / S
        qx = 0.25 * S
        qy = (R[0,1] + R[1,0]) / S
        qz = (R[0,2] + R[2,0]) / S
    elif (R[1,1] > R[2,2]):
        S = np.sqrt(1 + R[1,1] - R[0,0] - R[2,2]) * 2 # S=4*qy
        qw = (R[0,2] - R[2,0]) / S
        qx = (R[0,1] + R[1,0]) / S
        qy = 0.25 * S
        qz = (R[1,2] + R[2,1]) / S
    else:
        S = np.sqrt(1 + R[2,2] - R[0,0] - R[1,1]) * 2 # S=4*qz
        qw = (R[1,0] - R[0,1]) / S
        qx = (R[0,2] + R[2,0]) / S
        qy = (R[1,2] + R[2,1]) / S
        qz = 0.25 * S
    
    Q = [qw, qx, qy, qz]
    
    return Q

def spin_view(view,axis, angle, degrees=True):
    """特定の軸を中心に回転させる

    Parameters:
    view: nglview.NGLWidget
        NGLWidget
    axis: int
        0->x軸
        1->y軸
        2->z軸
    angle :float
        角度
    degrees: bool
        Tureの場合,Degree,Falseの場合Radian  
    """
    # NGLWidget.controlはViewerControlというクラスのapply_matrixメソッドを用いる
    if degrees:
        angle = radians(angle)
    if axis == 0:
        mat = [1,0,0,0,0,cos(angle),-sin(angle),0,0,sin(angle),cos(angle),0,0,0,1]
    elif axis == 1:
        mat = [cos(angle),0,sin(angle),0,0,1,0,0,-sin(angle),0,cos(angle),0,0,0,1]
    elif axis == 2:
        mat = [cos(angle),-sin(angle),0,0,sin(angle),cos(angle),0,0,0,0,1,0,0,0,1]
    else:
        raise("axis=1 or 2 or 3") 
    view.control.apply_matrix(mat)
