import sys
sys.path.insert(0, './')

import os
import json
import numpy as np
import open3d as o3d
import copy

import pathlib

current_path = str(pathlib.Path(__file__).parent.absolute())
path_mesh = current_path + '/obj/'
path_record = current_path + '/his/'

mesh_jsons = [f for f in os.listdir(path_mesh) if f.endswith('.' + 'json')]
mesh_names = [os.path.splitext(f)[0] for f in mesh_jsons]


record_jsons = [f for f in os.listdir(path_record) if f.endswith('.' + 'json')]
record_names = [os.path.splitext(f)[0] for f in record_jsons]
record_num = list(set([f.split('_')[0] for f in record_names]))
record_num = [int(f) for f in record_num]
record_num.sort()

out_mesh_dir = current_path + '/out1/'
os.makedirs(out_mesh_dir, exist_ok=True)
obj_dict = {}

for name in mesh_names:
    mesh_name = path_mesh + name + '.json'
    with open(mesh_name) as f:
        json_mesh = json.load(f)
    mesh_obj = json_mesh['object']['geometries'][0]['data']
    obj_name = out_mesh_dir + name + ".obj"
    print(obj_name)
    with open(obj_name, "w") as f:
        f.write(mesh_obj)
    mesh = o3d.io.read_triangle_mesh(obj_name)
    trans = np.array(json_mesh['object']['object']['matrix']).reshape((4,4)).transpose()
    mesh = copy.deepcopy(mesh).transform(trans)
    obj_dict[name] = mesh

for s in record_num:
    for name in mesh_names:
        record_name = path_record + str(s) + '_' + name + '.json'

        with open(record_name) as f:
            json_record = json.load(f)
        trans = np.array(json_record['matrix']).reshape((4,4)).transpose()
        mesh = obj_dict[name]
  
        t_mesh = copy.deepcopy(mesh).transform(trans)

        obj_name = out_mesh_dir + str(s) + "_" + name + ".obj"
        print(s, obj_name)
        o3d.io.write_triangle_mesh(obj_name, t_mesh)

