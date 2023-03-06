import sys
import traceback
from typing import Any,Dict,List,Optional,Union
import threading
import time
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io.vasp import write_vasp

from ase.io import Trajectory, write
from ase.visualize.ngl import NGLDisplay
from ipywidgets import (Button, Checkbox, Output,
                        Text, BoundedFloatText,RadioButtons,Image,ColorPicker,BoundedIntText,FloatSlider,
                        HBox,VBox,Tab,Dropdown,Layout)
from traitlets import Bunch

import nglview as nv
from nglview import NGLWidget
from nglview.component import ComponentViewer
from nglview.color import ColormakerRegistry

# USER
from grrmpy.visualize.functions import (update_tooltip_atoms,generate_js_code,
                                        get_struct,add_force_shape,rotate_view,spin_view)
from grrmpy.visualize.color import default,vesta,jmol

class View(NGLDisplay):
    """
    Parameters:
    
    atoms:
        | AtomsまたはAtomsのリストまたはTrajectoryクラス
    xsize: 
        | 横幅(px単位)
    ysize:
        | 縦幅(px単位)
    """
    def __init__(
        self,
        atoms: Union[Atoms, Trajectory, List[Atoms]],
        xsize: int = 400,
        ysize: int = 500,
        ):
        super().__init__(atoms, xsize=xsize, ysize=ysize)
        self.v = self.gui.view  # For backward compatibility...
        # del self.gui # デフォルトのGUIを削除
        # self.gui = HBox([self.view, VBox()]) # GUIを再設定
        
        # Make useful shortcuts for the user of the class
        self.gui.view = self.view
        self.gui.control_box = self.gui.children[1]
        self.gui.custom_colors = self.custom_colors
        
        ####### Property #####################################
        self.replace_structure = False
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
        #色の設定
        self.cm = ColormakerRegistry
        default_jscode = generate_js_code(default)
        vesta_jscode = generate_js_code(vesta)
        jmol_jscode = generate_js_code(jmol)
        self.cm.add_scheme_func('default',default_jscode)
        self.cm.add_scheme_func('vesta',vesta_jscode)
        self.cm.add_scheme_func('jmol',jmol_jscode)
        ###################################################

        # ---原子上にマウスを置いたときに,原子のindexと位置を表示する
        update_tooltip_atoms(self.view, self._get_current_atoms())
        
        # GUI作成&表示
        self.build_gui()
        
        # 初期表示
        self.view.camera = "orthographic" if self.camera_style=='平行投影' else "perspective"
        self.view.add_spacefill()
        self.view.add_ball_and_stick()
        self.view.add_label(
            color="black",
            labelType="text",
            labelText=["" for _ in range(len(self._get_current_atoms()))],
            zOffset=2.0,
            attachment="middle_center",
            radius=1,
        )
        self._update_repr()
        
        self.view.unobserve(NGLWidget._on_frame_changed)
        self.view.observe(self._on_frame_changed, names=["frame"])
        
    @property
    def camera_style(self):
        return self.gui.camera_radio_btn.value
            
    @property
    def show_charge(self):
        return self.gui.show_charge_checkbox.value
    
    @property
    def show_force(self):
        return self.gui.show_force_checkbox.value
    
    @property
    def transparent(self):
        return self.gui.transparent.value
    
    @property
    def factor(self):
        return self.gui.factor.value
    
    def build_gui(self):
        ###カメラ(RadioButton)###
        self.gui.camera_radio_btn = RadioButtons(
            options=['平行投影','透視投影'],
            value='平行投影',
            description='カメラ')
        self.gui.camera_radio_btn.observe(self.change_camera)
        ###ラベル(RadioButton)###
        self.gui.label_radio_btn = RadioButtons(
            options=['なし','インデックス','元素','電荷','FixAtoms'],
            value='なし',
            description='ラベル')
        self.gui.label_radio_btn.observe(self.change_label)
        ###セル(CheckBox)####
        self.gui.cell_check_box = Checkbox(value=True,description="セルユニット",)
        self.gui.cell_check_box.observe(self.show_unitcell)
        ###カラースキーム###
        self.csel = Dropdown(options=["default","vesta","element","jmol"],
                             value='element', description='色')
        self.csel.observe(self._update_repr)
        ###モデル###
        self.gui.model_radio_btn = RadioButtons(
            options=['球棒モデル','空間充填モデル'],
            value='球棒モデル',
            description='モデル')
        self.gui.model_radio_btn.observe(self._update_repr)
        ###再配置(チェックボックス)###
        self.gui.replace_structure_checkbox = Checkbox(
            value=self.replace_structure,
            description="再配置")
        self.gui.replace_structure_checkbox.observe(self.change_replace_structure)
        ###アウトプット(エラー表示)###
        self.gui.out_widget = Output(layout={"border": "0px solid black"})
        ##IMG##
        self.gui.a = Button(description='A',tooltip='A軸方向から見る',layout = Layout(width='30px'))
        self.gui.a.on_click(lambda e,x=180,y=-90,z=90: self.rotate_view(e,x,y,z))
        self.gui.b = Button(description='B',tooltip='B軸方向から見る',layout = Layout(width='30px'))
        self.gui.b.on_click(lambda e,x=90,y=0,z=-90: self.rotate_view(e,x,y,z))
        self.gui.c = Button(description='C',tooltip='C軸方向から見る',layout = Layout(width='30px'))
        self.gui.c.on_click(lambda e,x=180,y=0,z=0: self.rotate_view(e,x,y,z))
        self.gui.center = Button(tooltip='中心',layout = Layout(width='30px'),icon='bullseye')
        self.gui.center.on_click(lambda e,:self.center(e))
        # self.gui.up = Button(description='x',tooltip=' x軸周りで回転(上)',layout = Layout(width='30px'))
        # self.gui.up.on_click(lambda e,axis=0,angle=10: self.spin_view(e,axis,angle))
        # self.gui.down = Button(description='x*',tooltip='x軸周りで回転(下))',layout = Layout(width='30px'))
        # self.gui.left = Button(description='y*',tooltip='y軸周りで回転(左)',layout = Layout(width='30px'))
        # self.gui.right = Button(description='y',tooltip='y軸周りで回転(右)',layout = Layout(width='30px'))
        # self.gui.cc = Button(description='z',tooltip='z軸周りで回転(時計周り)',layout = Layout(width='30px'))
        # self.gui.rc = Button(description='z*',tooltip='z軸周りで回転(反時計周り)',layout = Layout(width='30px'))
        # self.gui.step = BoundedFloatText(value=10,min=0,max=360,step=10,
        #                                  description='Step(*):',
        #                                  layout = Layout(width='150px'))
        # ファイル名
        self.gui.transparent = Checkbox(value=True,description="PIG保存時に透過背景")
        if isinstance(self.atoms, Atoms):
            self.gui.filename_text = Text(value="Atoms", description="ファイル名: ",layout = Layout(width='200px'))
            self.extention_dict = {0:'traj', 1:'cif', 2:'xyz', 3:'html', 4:"vasp", 5:'png'}
            self.gui.file_extention = Dropdown(options=[(val,key) for key,val in self.extention_dict.items()],value=1,layout = Layout(width='70px'))
        else:# list of Atoms or Traj
            self.gui.filename_text = Text(value="Images", description="ファイル名: ",layout = Layout(width='200px'))
            self.extention_dict = {0:'traj',1:'traj+',2:'cif',3:'cif+',4:'xyz',5:'xyz+',6:'html',7:'html+',8:"vasp",9:'png'}
            self.gui.file_extention = Dropdown(options=[(val,key) for key,val in self.extention_dict.items()],value=1,layout = Layout(width='70px'))
        ##ダウンロード##
        self.gui.download = Button(description='PNGをダウンロード',
                                   tooltip='ローカルPCにPNGをダウンロードする')
        self.gui.download.on_click(self.download_image)
        ##保存##
        self.gui.save = Button(description='ファイルへ保存')
        self.gui.save.on_click(self.save_image)
        # 電荷表示
        self.gui.show_charge_checkbox = Checkbox(value=False,description="電荷:色",)
        self.gui.show_charge_checkbox.observe(self.show_charge_event)
        self.gui.charge_scale_slider = BoundedFloatText(
            value=1.0, min=0.0, max=100.0, step=0.1, description="電荷:スケール",style = {'description_width': 'initial'})
        self.gui.charge_scale_slider.observe(self.show_charge_event)
        self.gui.show_force_checkbox = Checkbox(value=False,description="力",)
        self.gui.show_force_checkbox.observe(self.show_force_event)
        self.gui.force_scale_slider = FloatSlider(
            value=0.5, min=0.0, max=100.0, step=0.1, description="力:スケール")
        self.gui.force_scale_slider.observe(self.show_force_event)
                
        # その他
        self.gui.color_picker = ColorPicker(concise=False,description='ラベルカラー:',value='black',
                                            layout = Layout(width='200px'),style = {'description_width': 'initial'})
        self.gui.color_picker.observe(self.update_label)
        self.gui.label_size = BoundedFloatText(value=1.0,min=0,max=3.0,step=0.1,description='ラベルサイズ:',
                                               layout = Layout(width='200px'),style = {'description_width': 'initial'})
        self.gui.label_size.observe(self.update_label)
        self.gui.charge_round = BoundedIntText(value=2,min=0,max=10,step=1,description='電荷 小数点以下:',
                                               layout = Layout(width='200px'),style = {'description_width': 'initial'})
        self.gui.charge_round.observe(self.change_label)
        self.gui.factor = BoundedIntText(value=4,min=0,max=30,step=1,description='PIGの解像度:',
                                               layout = Layout(width='200px'),style = {'description_width': 'initial'})
        ###表示####
        # r = list(self.gui.control_box.children)
        img1 = HBox([
            self.gui.a,
            self.gui.b,
            self.gui.c,
            self.gui.center,
            # self.gui.step　# 回転の中心を変える方法が分からないので非表示
            ])
        # img2 = HBox([
        #     self.gui.up,
        #     self.gui.down,
        #     self.gui.right,
        #     self.gui.left,
        #     self.gui.cc,
        #     self.gui.rc,
        #     ])
        general = VBox([
            img1,
            # img2, # 回転の中心を変える方法が分からないので非表示
            self.gui.transparent,
            HBox([self.gui.filename_text,self.gui.file_extention]),
            HBox([self.gui.download,self.gui.save]),
            self.gui.label_radio_btn,
            self.gui.show_charge_checkbox,
            self.gui.charge_scale_slider,
            self.gui.show_force_checkbox,
            self.gui.force_scale_slider,
            ])
        other = VBox([
            self.csel,
            self.rad,
            self.gui.model_radio_btn,
            self.gui.camera_radio_btn,
            self.gui.cell_check_box,
        ])
        detail = VBox([
            self.gui.color_picker,
            self.gui.label_size,
            self.gui.charge_round,
            self.gui.factor,
        ])
        
        self.tab = Tab([general,other,detail],_titles={0:"プロパティなど", 1:"スタイル",2:"その他"})
        self.gui.control_box.children = tuple([self.tab,
                                               self.gui.replace_structure_checkbox,
                                               self.gui.out_widget])
        
    def save_image(self,e=None):
        name = self.gui.filename_text.value
        atoms = self._get_current_atoms() # 現在表示中のAtoms
        if self.extention_dict[self.gui.file_extention.value] == "traj":
            write(f"{name}.traj",atoms,format="traj")
        elif self.extention_dict[self.gui.file_extention.value] == "traj+":
            write(f"{name}.traj",self.atoms,format="traj")
        elif self.extention_dict[self.gui.file_extention.value] == "cif":
            write(f"{name}.cif",atoms,format="cif")
        elif self.extention_dict[self.gui.file_extention.value] == "cif+":
            write(f"{name}.cif",self.atoms,format="cif")
        elif self.extention_dict[self.gui.file_extention.value] == "xyz":
            write(f"{name}.xyz",atoms)
        elif self.extention_dict[self.gui.file_extention.value] == "xyz+":
            write(f"{name}.xyz",self.atoms)
        elif self.extention_dict[self.gui.file_extention.value] == "html":
            nv.write_html(f"{name}.html",self.view,tuple([self.view.frame]))
        elif self.extention_dict[self.gui.file_extention.value] == "html+":
            nv.write_html(f"{name}.html",self.view,(0,len(self.atoms)-1))
        elif self.extention_dict[self.gui.file_extention.value] == "vasp":
            write_vasp(f"{name}",atoms,sort=True,wrap=True,direct=True)
        elif self.extention_dict[self.gui.file_extention.value] == "png":
            thread = threading.Thread(target=self._save_image_png, args=(f"{name}.png", self.view), daemon=True)
            thread.start()
        
    def download_image(self,e=None):
        try:
            filename = self.gui.filename_text.value
            self.view.download_image(filename=filename,transparent=self.transparent,factor=self.factor)
        except Exception as e:
            with self.gui.out_widget:
                print(traceback.format_exc(), file=sys.stderr)

    def _save_image_png(self,filename: str, v: NGLWidget):
        try:
            image = v.render_image(transparent=self.transparent,factor=self.factor)
        except Exception as e:
            with self.gui.out_widget:
                print(traceback.format_exc(), file=sys.stderr)
        while not image.value:
            time.sleep(0.1)
        with open(filename, "wb") as fh:
            fh.write(image.value)
        
    def rotate_view(self,e,x,y,z):
        rotate_view(self.view,x=x,y=y,z=z)
        self.view.center()
        
    def spin_view(self,e, axis, angle):
        spin_view(self.view, axis, angle)
        self.view.center()
        
    def center(self,e=None):
        self.view.center()
        
    def update_label(self,e=None):
        self.view.update_label(
            color=self.gui.color_picker.value,
            labelType="text",
            labelText=self.labelText,
            zOffset=2.0,
            attachment="middle_center",
            radius=self.gui.label_size.value,
        )
        
    def _change_label(self,atoms,option):
        if option == "インデックス":
            self.labelText=[str(i) for i in range(len(atoms))]
        elif option == "元素":
            self.labelText=[i for i in atoms.get_chemical_symbols()]
        elif option == "電荷":
            try:
                self.labelText = np.round(atoms.get_charges().ravel(), self.gui.charge_round.value).astype("str").tolist()
            except:
                self.labelText=["" for _ in range(len(atoms))]
                with self.gui.out_widget:
                    raise Exception("Calculatorを設定してください")    
        elif option == "FixAtoms":
            self.labelText = self._get_fix_atoms_label_text(atoms)
        self.update_label()
        
    def change_label(self,e=None):
        self.gui.out_widget.clear_output()
        option = self.gui.label_radio_btn.value
        atoms = self._get_current_atoms()
        if option == "なし":
            self.labelText=["" for _ in range(len(atoms))]
            self.update_label()
            return 
        else:
            self._change_label(atoms,option)
            return 
           
    def change_camera(self,e=None):
        option = self.gui.camera_radio_btn.value
        if option == "平行投影":
            self.view.camera = 'orthographic'
        elif option == "透視投影":
            self.view.camera = 'perspective'
            
    def _update_repr(self,e=None):
        option = self.gui.model_radio_btn.value
        if option == "球棒モデル":
            self.view.update_spacefill(radiusType='covalent',
                                    radiusScale=self.rad.value,
                                    color_scheme=self.csel.value)#color_scale='rainbow')
            self.view.update_ball_and_stick(color_scheme=self.csel.value)
        elif option == "空間充填モデル":
            self.view.remove_spacefill()
            self.view.add_spacefill()
            self.view.update_spacefill(radiusType="vwf",color_scheme=self.csel.value)
                
    def show_unitcell(self,e=None):
        if self.gui.cell_check_box.value: 
            self.view.add_unitcell() # Cellの表示
        else:
            self.view.remove_unitcell()
                
    def _get_current_atoms(self) -> Atoms:
        if isinstance(self.atoms, Atoms):
            return self.atoms
        else:
            return self.atoms[self.view.frame]
        
    def _get_fix_atoms_label_text(self,atoms):
        indices_list = []
        for constraint in atoms.constraints:
            if isinstance(constraint, FixAtoms):
                indices_list.extend(constraint.index.tolist())
        label_text = []
        for i in range(len(atoms)):
            if i in indices_list:
                label_text.append("Fix")
            else:
                label_text.append("")
        return label_text
        
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
        
        ###ラベルの更新
        option = self.gui.label_radio_btn.value
        if option != "なし":
            atoms = self._get_current_atoms()
            self._change_label(atoms,option)
            
        
    def _ipython_display_(self, **kwargs):
        """viewプロパティを書かなくてもjupyter上で勝手に表示してくれる"""
        return self.gui._ipython_display_(**kwargs)
        
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
            
    def show_force_event(self, event: Optional[Bunch] = None):
        self.gui.out_widget.clear_output()
        self.clear_force()
        if self.show_force:
            self.add_force()
            
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
        
    def clear_force(self):
        # Remove existing force components.
        for c in self._force_components:
            self.v.remove_component(c)  # Same with c.clear()
        self._force_components = []
        
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

