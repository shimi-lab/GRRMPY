import sys
import traceback
from typing import Any, Dict, List, Optional, Union
from io import StringIO
import threading
import time

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io.proteindatabank import write_proteindatabank
from ase.io import Trajectory, write
from ase.visualize.ngl import NGLDisplay
from ipywidgets import Button, Checkbox, HBox, Output, Text, BoundedFloatText,RadioButtons
from traitlets import Bunch

import nglview as nv
from nglview import NGLWidget
from nglview.component import ComponentViewer


class View2(NGLDisplay):
    """
        View2の機能は全てgrrmpy.visualize.custum_viewer.Viewに移行済み
        このクラスは将来的に消す予定
    
        Parameters:
        
        atoms: Union[Atoms, Trajectory, List[Atoms]]
            Atomsオブジェクト,TrajectoryまたはAtomsリスト
        replace_structure: bool
            再配置する
        representations: list
            デフォルトは["ball+stick"]
        xsize: integer
            横幅
        ysize: int
            縦幅
        show_force: bool
            力のチェックボックスの初期値
        show_charge: bool
            電荷のチェックボックスの初期値
            
        Example:
        
        >>> from grrmpy.visualize import View
        >>> View(atoms)
    
    """

    def __init__(
        self,
        atoms: Union[Atoms, Trajectory, List[Atoms]],
        replace_structure: bool = False,
        representations = ["ball+stick"],
        xsize: int = 500,
        ysize: int = 500,
        show_force: bool = False,
        show_charge: bool = False,
    ):
        super().__init__(atoms, xsize=xsize, ysize=ysize)
        self.v = self.gui.view  # For backward compatibility...
        self.view = self.gui.view
        for rep in representations:
            self.view.add_representation(rep)
        # --- Attribute ---
        self.replace_structure = replace_structure
        self._use_struct_cache = True
        self._struct_cache = []
        self._force_components = []
        self.force_color = [1, 0, 0]  # Red vector for force color.
        self.pre_label = False

        if isinstance(atoms, Atoms):
            self._struct_cache = [None]
        else:
            # atoms is Trajectory or List[Atoms]
            self._struct_cache = [None for _ in range(len(atoms))]

        # --- NGL Viewer custom setting ---
        update_tooltip_atoms(self.view, self._get_current_atoms())

        # --- GUI ---
        self.gui.save_current_atoms_checkbox = Checkbox(
            value=False,
            description="現在の構造のみ保存",
        )
        self.gui.filename_text = Text(value="sample.traj", description="ファイル名: ")
        self.gui.download_image_button = Button(
            description="Pngをダウンロード",
            tooltip="現在表示中の構造をPNGとしてローカルPCに直接保存する",
        )
        self.gui.download_image_button.on_click(self.download_image)
        self.gui.save_image_button = Button(
            description="保存",
            tooltip="ファイルへ保存する.\n"
            ".png,.htmlとその他ASEのwriteで指定できる拡張子で保存できる.\n"
        )
        self.gui.save_image_button.on_click(self.save_image)
        # Show force vector. Note that Calculator must be set to atoms.
        self.gui.show_force_checkbox = Checkbox(
            value=False,
            description="力",
        )
        self.gui.show_force_checkbox.observe(self.show_force_event)
        self.gui.force_scale_slider = BoundedFloatText(
            value=0.5, min=0.0, max=100.0, step=0.1, description="力:スケール"
        )
        self.gui.force_scale_slider.observe(self.show_force_event)

        # Show charge color. Note that Calculator must be set to atoms.
        self.gui.show_charge_checkbox = Checkbox(
            value=False,
            description="電荷:色",
        )
        self.gui.show_charge_checkbox.observe(self.show_charge_event)
        self.gui.charge_scale_slider = BoundedFloatText(
            value=1.0, min=0.0, max=100.0, step=0.1, description="電荷:スケール"
        )
        self.gui.charge_scale_slider.observe(self.show_charge_event)
         
        self.gui.label_radiobutton = RadioButtons(
            options=['なし', '電荷', 'インデックス'],
            value='なし', # Defaults to 'pineapple'
            description='ラベル',
            disabled=False
            )
        self.gui.label_radiobutton.observe(self.change_label)
        
        # Now the border is 0px and it is invisible at first. Text is visible only when text is written.
        self.gui.out_widget = Output(layout={"border": "0px solid black"})
        self.gui.replace_structure_checkbox = Checkbox(
            value=self.replace_structure,
            description="再配置",
        )
        self.gui.replace_structure_checkbox.observe(self.change_replace_structure)
        

        r = list(self.gui.control_box.children)
        r += [
            self.gui.save_current_atoms_checkbox,
            self.gui.filename_text,
            HBox(
                [
                    self.gui.download_image_button,
                    self.gui.save_image_button,
                ]
            ),
            self.gui.show_force_checkbox,
            self.gui.force_scale_slider,
            self.gui.show_charge_checkbox,
            self.gui.charge_scale_slider,
            self.gui.label_radiobutton,
            self.gui.replace_structure_checkbox,
            self.gui.out_widget,
        ]
        self.gui.control_box.children = tuple(r)

        # Unsubscribe original `_on_frame_changed`, and subscribe custom `_on_frame_changed`.
        self.view.unobserve(NGLWidget._on_frame_changed)
        self.view.observe(self._on_frame_changed, names=["frame"])

        # Set & apply show_force, show_charge effect
        self.gui.show_force_checkbox.value = show_force
        self.gui.show_charge_checkbox.value = show_charge

    def change_label(self, event: Optional[Bunch] = None):
        self.gui.out_widget.clear_output()
        option = self.gui.label_radiobutton.value
        if option != "なし" and not self.pre_label:
            self.show_label(option)
        elif option != "なし" and self.pre_label:
            self.update_label(option)
        else:
            self.view.remove_label()
        self.pre_label = True if option!="なし" else False

    def download_image(self, clicked_button: Optional[Button] = None):
        filename = self.gui.filename_text.value
        self.view.download_image(filename=filename)

    def save_image(self, clicked_button: Optional[Button] = None):
        filename = self.gui.filename_text.value
        format = "vasp" if filename=="POSCAR" or filename=="CONTCAR" else None
        if self.save_current_atoms:
            if any([filename.endswith(".png"),filename.endswith(".html")]):
                save_image(filename, self.view)
            else:
                write(filename,self._get_current_atoms(),format=format)
        else:
            if filename.endswith(".png"):
                f"PNGの場合,現在の構造のみを保存するにチャックを入れて下さい．"
            elif filename.endswith(".html"):
                nv.write_html(filename,[self.view],(0,len(self.atoms)-1))
            else:
                write(filename,self.atoms,format=format)

    def _get_current_atoms(self) -> Atoms:
        if isinstance(self.atoms, Atoms):
            return self.atoms
        else:
            return self.atoms[self.view.frame]

    def show_charge_event(self, event: Optional[Bunch] = None, refresh: bool = True):
        self.gui.out_widget.clear_output()
        if self.show_charge:
            atoms = self._get_current_atoms()
            # TODO: How to change `scale` and `radiusScale` by user?
            # Register "atomStore.partialCharge" attribute inside javascript
            charge_scale: float = self.gui.charge_scale_slider.value
            # Note that Calculator must be set here!
            try:
                charge_str = str((atoms.get_charges().ravel() * charge_scale).tolist())
            except Exception as e:
                with self.gui.out_widget:
                    print(traceback.format_exc(), file=sys.stderr)
                # `append_stderr` method shows same text twice somehow...
                # self.gui.out_widget.append_stderr(str(e))
                return

            var_code = f"var chargeArray = {charge_str}"
            js_code = """
            var component = this.stage.compList[0]
            var atomStore = component.structure.atomStore
            if (atomStore.partialCharge === undefined) {
                atomStore.addField('partialCharge', 1, 'float32')
            }

            for (let i = 0; i < chargeArray.length; ++i) {
              atomStore.partialCharge[i] = chargeArray[i];
            }
            """
            self.view._execute_js_code(var_code + js_code)

            # Show charge color
            # TODO: More efficient way:
            #  We must set other "color_scheme" at first, to update "partialcharge" color scheme...
            # color_schme="element" is chosen here, but any color_scheme except "partialcharge" is ok.
            # Skip this procedure to avoid heavy process, user must turn on and off "show charge" now.
            if refresh:
                self.view._update_representations_by_name(
                    "spacefill",
                    radiusType="covalent",
                    radiusScale=self.rad.value,
                    color_scheme="element",
                    color_scale="rwb",
                )
            self.view._update_representations_by_name(
                "spacefill",
                radiusType="covalent",
                radiusScale=self.rad.value,
                color_scheme="partialcharge",
                color_scale="rwb",
            )
        else:
            # Revert to original color scheme.
            self._update_repr()
            
    def change_replace_structure(self,event: Optional[Bunch] = None):
        if self.gui.replace_structure_checkbox.value:
            self.replace_structure = True
            self._on_frame_changed(None)
        else:
            self.replace_structure = False

    def _on_frame_changed(self, change: Dict[str, Any]):
        """set and send coordinates at current frame"""
        v: NGLWidget = self.view
        atoms: Atoms = self._get_current_atoms()

        self.clear_force()
        if self.replace_structure:
            # set and send coordinates at current frame
            struct = self._struct_cache[v.frame]
            if struct is None:
                struct = get_struct(atoms)
                if self._use_struct_cache:
                    self._struct_cache[v.frame] = struct  # Cache
            v._remote_call("replaceStructure", target="Widget", args=struct)
        else:
            # Only update position info
            v._set_coordinates(v.frame)

        if self.show_force:
            self.add_force()
        if self.show_charge:
            self.show_charge_event()

        # Tooltip: update `var atoms_pos` inside javascript.
        atoms = self._get_current_atoms()
        if atoms.get_pbc().any():
            _, Q = atoms.cell.standard_form()
        else:
            Q = np.eye(3)
        Q_str = str(Q.tolist())
        var_str = f"this._Q = {Q_str}"
        v._execute_js_code(var_str)

    def clear_force(self):
        # Remove existing force components.
        for c in self._force_components:
            self.v.remove_component(c)  # Same with c.clear()
        self._force_components = []

    def add_force(self):
        force_scale: float = self.gui.force_scale_slider.value
        try:
            atoms = self._get_current_atoms()
            c = add_force_shape(atoms, self.v, force_scale, self.force_color)
            self._force_components.append(c)
        except Exception as e:
            with self.gui.out_widget:
                print(traceback.format_exc(), file=sys.stderr)
            # `append_stderr` method shows same text twice somehow...
            # self.gui.out_widget.append_stderr(str(e))
            return

    def show_force_event(self, event: Optional[Bunch] = None):
        self.gui.out_widget.clear_output()
        self.clear_force()
        if self.show_force:
            self.add_force()

    @property
    def show_force(self):
        return self.gui.show_force_checkbox.value

    @property
    def show_charge(self):
        return self.gui.show_charge_checkbox.value
        
    @property
    def save_current_atoms(self):
        return self.gui.save_current_atoms_checkbox.value
    
    def _ipython_display_(self, **kwargs):
        """Called when `IPython.display.display` is called on the widget."""
        return self.gui._ipython_display_(**kwargs)

    def show_label(self, option, color: str = "black", radius: float = 0.8):
        atoms = self._get_current_atoms()
        if option == "電荷":
            try:
                labelText = np.round(atoms.get_charges().ravel(), 3).astype("str").tolist()
            except Exception as e:
                with self.gui.out_widget:
                    print(traceback.format_exc(), file=sys.stderr)
                # `append_stderr` method shows same text twice somehow...
                # self.gui.out_widget.append_stderr(str(e))
                return 
        elif option == "インデックス":
            labelText = [str(i) for i in range(len(atoms))]
        self.view.add_label(
            color=color,
            labelType="text",
            labelText=labelText,
            zOffset=2.0,
            attachment="middle_center",
            radius=radius,
        )
        
    def update_label(self, option, color: str = "black", radius: float = 0.8):
        """option:`電荷`or`インデックス`"""
        atoms = self._get_current_atoms()
        if option == "電荷":
            try:
                labelText = np.round(atoms.get_charges().ravel(), 3).astype("str").tolist()
            except Exception as e:
                with self.gui.out_widget:
                    print(traceback.format_exc(), file=sys.stderr)
                # `append_stderr` method shows same text twice somehow...
                # self.gui.out_widget.append_stderr(str(e))
                return 
        elif option == "インデックス":
            labelText = [str(i) for i in range(len(atoms))]
        self.view.update_label(
            color=color,
            labelType="text",
            labelText=labelText,
            zOffset=2.0,
            attachment="middle_center",
            radius=radius,
        )

    def _get_fix_atoms_label_text(self):
        atoms = self._get_current_atoms()
        indices_list = []
        for constraint in atoms.constraints:
            if isinstance(constraint, FixAtoms):
                indices_list.extend(constraint.index.tolist())
        label_text = []
        for i in range(len(atoms)):
            if i in indices_list:
                label_text.append("FixAtoms")
            else:
                label_text.append("")
        return label_text

    def show_fix_atoms_constraint(self, color: str = "black", radius: float = 0.8):
        self.view.add_label(
            color=color,
            labelType="text",
            labelText=self._get_fix_atoms_label_text(),
            zOffset=2.0,
            attachment="middle_center",
            radius=radius,
        )

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

def _save_image_png(filename: str, v: NGLWidget):
    """Save nglview image in png format.

    Note that it should be run on another thread.
    See: https://github.com/nglviewer/nglview/blob/master/docs/FAQ.md#how-to-make-nglview-view-object-write-png-file

    Args:
        filename (str):
        v (NGLWidget):
    """
    image = v.render_image()
    while not image.value:
        time.sleep(0.1)
    with open(filename, "wb") as fh:
        fh.write(image.value)


def save_image(filename: str, v: NGLWidget):
    """Save nglview image based on file name extension.

    Args:
        filename (str):
        v (NGLWidget):
    """        
    if filename.endswith(".png"):
        thread = threading.Thread(
            target=_save_image_png, args=(filename, v), daemon=True
        )
        thread.start()
    elif filename.endswith(".html"):
        nv.write_html(filename, [v])  # type: ignore
    else:
        print(f"filename {filename}: extension not supported!")
        
def update_tooltip_atoms(view: NGLWidget, atoms: Atoms):
    """Update Tooltip text to show atom index & positions when mousehover on atoms

    Args:
        view:
        atoms:

    Returns:

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
    # Above code overrides `stage.tooltip`.
    # It is used in `_onMove` method to remove tooltip when un-hovered.
    # https://github.com/nglviewer/ngl/blob/bd4a31c72e007d170b6bae298a5f7c976070e173/src/stage/mouse-behavior.ts#L31-L33
    view._execute_js_code(var_str + script_str)