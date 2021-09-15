#!/usr/bin/python
#
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Lint as: python3
"""Example how to use differentiable physics engine using python.
We used pybind11 to expose the classes of pydiffphys.
At the moment, only double precision version is exposed.
Will also expose stan_math forward mode differentiation, dual and fixed point.
"""

from absl import app
from absl import flags

import torch
import pydiffarti as pd
import api_diff
import numpy as np

FLAGS = flags.FLAGS



def get_loss(joints, target_joint, param_ctl, q):
  n_link = 1
  n_start = -3*n_link
  if n_link == 1:
    ours = joints[n_start:]
  else:
    ours = joints[n_start:n_start+3]
  target = target_joint[-3:]
  loss = torch.sqrt(((ours[-3:]-target[-3:])**2).mean())
  loss_reg = torch.norm(param_ctl) * 0.1
  return loss #+ loss_reg 



def construct_ini_tau(param_ctl, n_step):
  hold_ratio = 0.6
  n_hold = int(hold_ratio * n_step)
  return param_ctl



def main(argv):
  if len(argv) > 1:
    raise app.UsageError("Too many command-line arguments.")

  world = pd.TinyWorld(do_vis=True)
  world.friction = pd.Utils.scalar_from_double(2)
  world.restitution = pd.Utils.scalar_from_double(0.5)
  parser = pd.TinyUrdfParser()
  convert_tool = pd.UrdfToMultiBody2()


  mb_gripper = world.create_multi_body()
  mb_gripper.isFloating = False
  urdf_structures = parser.load_urdf('./data/punch/p2.urdf')
  convert_tool.convert2(urdf_structures, world, mb_gripper)
  ini_pos = pd.TinyVector3(pd.Utils.scalar_from_double(0), pd.Utils.scalar_from_double(0), pd.Utils.scalar_from_double(0))
  mb_gripper.set_position(ini_pos)
  print(mb_gripper.dof_u())
  print(mb_gripper.dof())
  print(mb_gripper.dof_qd())
  print(mb_gripper.dof_state())

  world.sync_visual_meshcat(0)

  mb_sphere = world.create_multi_body()
  mb_sphere.isFloating = True
  urdf_structures = parser.load_urdf('./data/punch/sphere.urdf')
  convert_tool.convert2(urdf_structures, world, mb_sphere)

  world.sync_visual_meshcat(0)


  knee_angle = -0.5
  abduction_angle = 0.2
  init_q = torch.tensor([0,0,
                         1,0,0,0,0,-0.1,-0.85], dtype=torch.float32, requires_grad=True)
  init_qd = torch.tensor([0,0,
                          0,0,0,0,0,0], dtype=torch.float32, requires_grad=False)

  grav = pd.TinyVector3(pd.Utils.zero(), pd.Utils.zero(), pd.Utils.fraction(-980, 100))

  dt = 1/1000
  finer = 1
  world.dt = pd.Utils.scalar_from_double(dt/finer)
  n_step = 1000
  dof_u = 2
  n = dof_u * n_step
  dim_control = 2

  param_ctl = torch.normal(mean=0, std=1,
    size=(n_step, dim_control), dtype=torch.float32, requires_grad=True)

  target_joint = [ 0.0000, -3, -1]
  target_joint = torch.tensor(target_joint, dtype=torch.float32, requires_grad=True)

  world.adj_initialize(grav, n_step, dof_u)
  optimizer = torch.optim.Adam([param_ctl], lr=0.1)
  frameskip_gfx_sync = 16

  for steps in range(1000):
    init_tau = construct_ini_tau(param_ctl, n_step)

    sync_counter = frameskip_gfx_sync
    optimizer.zero_grad()
    q, qd = init_q, init_qd
    for sim_step in range(n_step):
      tau = init_tau[sim_step]

      for i in range(finer):
        q, qd = api_diff.sim_layer(q, qd, tau, world)

      if steps % 10 == 0:
        sync_counter = sync_counter+1
        if sync_counter > frameskip_gfx_sync:
          world.sync_visual_meshcat(sim_step+1)
          sync_counter = 0

    joints = api_diff.get_joints(q, world)
    
    loss = get_loss(joints, target_joint, param_ctl, q)
    print(steps, " loss =",loss)
    if loss.detach().numpy() < 5e-2:
      print(steps)
      break
    loss.backward()
    optimizer.step()
  print("done")
  # print('tau=',init_tau.reshape([-1,5]))
  # print(init_tau.min(), init_tau.max())
  # np.save('./data/push/tau.npy', init_tau.reshape([-1,5]).detach().numpy())
  while True:
    continue


if __name__ == "__main__":
  app.run(main)
