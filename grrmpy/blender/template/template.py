import bpy
import numpy as np
from mathutils import Matrix
import zipfile
import pickle
import json
from pathlib import Path


data_list = {{data_list}}


def delete_all_objects():
    for col in bpy.data.collections:
        for item in col.objects:
            col.objects.unlink(item)
            bpy.data.objects.remove(item)
    for item in bpy.context.scene.collection.objects:
        bpy.context.scene.collection.objects.unlink(item)
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)
    for item in bpy.data.meshes:
        bpy.data.meshes.remove(item)
    for item in bpy.data.materials:
        bpy.data.materials.remove(item)
    return

{%- for data in data_list %}
{%- if data.get("subdivision_surface",False) %}

def apply_subdivision_surface(obj):
    obj.modifiers.new("subd", type='SUBSURF')
    obj.modifiers['subd'].levels = {{data["subdivision_surface"]["level"]}}
    obj.modifiers['subd'].render_levels = {{data["subdivision_surface"]["render_levels"]}}
{% break %}
{%- endif %}
{%- endfor %}

{%- for data in data_list %}
{%- if not data["style"] in ["space_filling","animation"] %}
def distance(a, b):
    return np.sqrt(np.dot(a - b, a - b))
 
def normalize(a):
    return np.array(a) / np.sqrt(np.dot(a, a))

def rotate_object(obj,rot_mat):
    orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
    orig_loc_mat   = Matrix.Translation(orig_loc)
    orig_rot_mat   = orig_rot.to_matrix().to_4x4()
    orig_scale_mat = (Matrix.Scale(orig_scale[0],4,(1,0,0)) @ 
                      Matrix.Scale(orig_scale[1],4,(0,1,0)) @ 
                      Matrix.Scale(orig_scale[2],4,(0,0,1)))
    obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat 
{% break %}
{%- endif %}
{%- endfor %}

{%- for data in data_list %}
{%- if not data["style"] in ["stick"] %}
def draw_atoms(name, elements, positions, ball_sizes, subdivision_surface):
    for i,(element, position) in enumerate(zip(elements, positions)):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=ball_sizes[element], location=position)
        bpy.context.active_object.data.materials.append(bpy.data.materials[f"{name}{element}"])
        bpy.context.active_object.name = f"{name}Atom{i}{element}"
        bpy.ops.object.shade_smooth()
        if subdivision_surface:
            apply_subdivision_surface(bpy.context.active_object)
{% break %}
{%- endif %}
{%- endfor %}

{%- for data in data_list %}
{%- if data["style"] in ["stick","ball_and_stick"]%}
{%- if not data.get("bicolor",False)%}
def draw_mono_color_bonds(name,bonds,positions,bond_radius):
    for atom_1, atom_2 in bonds:
        pos_1 = positions[atom_1]
        pos_2 = positions[atom_2]
            
        difference = pos_2 - pos_1
        center = (pos_2 + pos_1) / 2.0
        magnitude = distance(pos_1, pos_2)
            
        bond_direction = normalize(difference)
        vertical = np.array((0.0, 0.0, 1.0))
        rotation_axis = np.cross(bond_direction, vertical)
        angle = -np.arccos(np.dot(bond_direction, vertical))
    
        bpy.ops.mesh.primitive_cylinder_add(radius=bond_radius, 
                                            depth=magnitude, 
                                            location=center)
        bpy.context.active_object.data.materials.append(bpy.data.materials[f'{name}bond'])
        bpy.context.active_object.name = f"{name}Bond({atom_1}-{atom_2}){i}"
        bpy.ops.object.shade_smooth()
        rotate_object(bpy.context.active_object, Matrix.Rotation(angle, 4, rotation_axis))
{% break %}
{%- endif %}
{%- endif %}
{%- endfor %}

{%- for data in data_list %}
{%- if data["style"] in ["stick","ball_and_stick"]%}
{%- if data.get("bicolor",False)%}
def draw_bicolor_bonds(bonds,positions,elements,bond_radius,half):
    for atom_1, atom_2 in bonds:
        pos_1 = positions[atom_1]
        pos_2 = positions[atom_2]
        elements_1 = elements[atom_1]
        elements_2 = elements[atom_2]
        
        if half:
            difference = pos_2 - pos_1
            position_1 = (pos_2*1 + pos_1*3) / 4
            position_2 = (pos_2*3 + pos_1*1) / 4
            magnitude1 = distance(pos_1, pos_2)/2
            magnitude2 = magnitude1
        else:
            size_1 = ball_sizes[elements_1]
            size_2 =  ball_sizes[elements_2]
            difference = pos_2 - pos_1
            d = distance(pos_1, pos_2)
            l = (d- size_1 - size_2)/2
            ratio_1 = (size_1+l)/d
            ratio_2 = (size_2+l)/d
            position_1 = (pos_2*(ratio_1*0.5) + pos_1*(1-(ratio_1*0.5)))
            position_2 = (pos_2*(1-(ratio_2*0.5)) + pos_1*(ratio_2*0.5))
            magnitude1 = d*ratio_1
            magnitude2 = d*ratio_2

                
        bond_direction = normalize(difference)
        vertical = np.array((0.0, 0.0, 1.0))
        rotation_axis = np.cross(bond_direction, vertical)
        angle = -np.arccos(np.dot(bond_direction, vertical))
    
        bpy.ops.mesh.primitive_cylinder_add(radius=bond_radius, 
                                            depth=magnitude1, 
                                            location=position_1)
        bpy.context.active_object.data.materials.append(bpy.data.materials[f"{name}{elements_1}"])
        bpy.context.active_object.name = f"{name}Bond({atom_1}-{atom_2}){i}"
        bpy.ops.object.shade_smooth()
        rotate_object(bpy.context.active_object, Matrix.Rotation(angle, 4, rotation_axis))
        
        bpy.ops.mesh.primitive_cylinder_add(radius=bond_radius, 
                                            depth=magnitude2, 
                                            location=position_2)                                            
        bpy.context.active_object.data.materials.append(bpy.data.materials[f"{name}{elements_2}"])
        bpy.context.active_object.name = f"{name}Bond({atom_2}-{atom_1}){i}"
        bpy.ops.object.shade_smooth()
        rotate_object(bpy.context.active_object, Matrix.Rotation(angle, 4, rotation_axis))
{% break %}
{%- endif %}
{%- endif %}
{%- endfor %}

{%- for data in data_list %}
{%- if data["style"] =="animation" %}
def add_keyflame(name,frame_num,step,positions,chemical_symbols):
    bpy.context.scene.frame_set(frame_num)
    for i,(position,element) in enumerate(zip(positions,chemical_symbols)):
        obj = bpy.data.objects[f"{name}Atom{i}{element}"]
        obj.location = position
        obj.keyframe_insert(data_path = 'location',index = -1)
    frame_num += step
    return frame_num
{% break %}
{%- endif %}
{%- endfor %}

def register_materials(name,rgba,cartoon):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    if cartoon["apply"]:
        snf = mat.node_tree.nodes.new("ShaderNodeFresnel")
        snmr = mat.node_tree.nodes.new("ShaderNodeMixRGB")
        bsdf.inputs[7].default_value = 1
        snmr.inputs['Color1'].default_value = rgba
        snmr.inputs['Color2'].default_value = cartoon["color"]
        snf.inputs['IOR'].default_value = cartoon["IOR"]
        snf.location = -500, -110
        snmr.location = -300, -110
        mat.node_tree.links.new(snf.outputs[0], snmr.inputs[0])
        mat.node_tree.links.new(snmr.outputs[0], bsdf.inputs[0])
    else:
        bsdf.inputs[0].default_value = rgba

delete_all_objects()
for i,data in enumerate(data_list):
    name = "" if len(data_list)==1 else f"{i}_"
    for symb,rgba in data["colors"].items():
        register_materials(f"{name}{symb}",rgba,cartoon=data["cartoon"])
        
    if data["style"] == "animation":
        ball_sizes = {symb:data["scale"]*size for symb,size in data["sizes"].items()}
        step = data["step"]
        start = data["start"]
        frame_num = start
        p = Path(bpy.data.filepath)
        pkl_path = str(p.with_name(data["file"]).resolve())
        with open(pkl_path,"rb") as f:
            positions = pickle.load(f)
            ball_sizes = {symb:data["scale"]*size for symb,size in data["sizes"].items()}
            subdivision_surface = data["subdivision_surface"]["apply"]
            draw_atoms(name,data["chemical_symbols"],np.array(positions),ball_sizes,subdivision_surface)
            frame_num = add_keyflame(name,frame_num,step,positions,data["chemical_symbols"])
            while True:
                try:
                    positions = pickle.load(f)
                    frame_num = add_keyflame(name,frame_num,step,positions,data["chemical_symbols"])
                except EOFError:
                    break
        continue
    
    if "stick_color" in data.keys():
        register_materials(f"{name}bond",rgba=data["stick_color"],cartoon=data["cartoon"])
    
    positions = np.array(data["positions"])
    if data["style"] in ["ball_and_stick","space_filling","animation"]:
        ball_sizes = {symb:data["scale"]*size for symb,size in data["sizes"].items()}
        subdivision_surface = data["subdivision_surface"]["apply"]
        draw_atoms(name,data["chemical_symbols"],positions,ball_sizes,subdivision_surface)
    if data["style"] in ["stick","ball_and_stick"]:
        if data["bicolor"]:
            half = True if data["style"] == "stick" else False
            draw_bicolor_bonds(data["bonds"],positions,data["chemical_symbols"],data["radius"],half=half)
        else:
            draw_mono_color_bonds(name,data["bonds"],positions,data["radius"])

        

        
    


