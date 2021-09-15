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

def get_loss(joints, target_joint, param_ctl):
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
  hold_ratio = 0.9
  n_hold = int(hold_ratio * n_step)
  # n_hold = n_step


  init_tau = []
  const_input = torch.tensor([-3e-3,-3e-3], requires_grad=False) * 0
  const_hold = const_input.repeat(n_hold, 1)
  const_release = const_input.repeat(n_step - n_hold, 1) * 0
  const_input = torch.cat([const_hold, const_release], 0)

  init_tau = torch.cat([param_ctl, const_input], dim=1)
 
  return init_tau

def main(argv):
  if len(argv) > 1:
    raise app.UsageError("Too many command-line arguments.")

  world = pd.TinyWorld(do_vis=True)
  print(pd.Utils.getDouble(world.friction))
  world.friction = pd.Utils.scalar_from_double(2)
  print(pd.Utils.getDouble(world.friction))
  parser = pd.TinyUrdfParser()
  convert_tool = pd.UrdfToMultiBody2()


  mb_gripper = world.create_multi_body()
  mb_gripper.isFloating = False
  urdf_structures = parser.load_urdf('./data/franka_panda/panda.urdf')
  convert_tool.convert2(urdf_structures, world, mb_gripper)
  ini_pos = pd.TinyVector3(pd.Utils.scalar_from_double(0), pd.Utils.scalar_from_double(0), pd.Utils.scalar_from_double(0))
  mb_gripper.set_position(ini_pos)
  print(mb_gripper.dof_u())
  print(mb_gripper.dof())
  print(mb_gripper.dof_qd())
  print(mb_gripper.dof_state())


  mb_sphere = world.create_multi_body()
  mb_sphere.isFloating = True
  urdf_structures = parser.load_urdf('./data/franka_panda/sphere2.urdf')
  convert_tool.convert2(urdf_structures, world, mb_sphere)
  ini_pos = pd.TinyVector3(
    pd.Utils.scalar_from_double(9.0644e-02), 
    pd.Utils.scalar_from_double(-1.1667e-02), 
    pd.Utils.scalar_from_double(8.5423e-01))
  mb_sphere.set_position(ini_pos)



  knee_angle = -0.5
  abduction_angle = 0.2

  init_q = torch.tensor([0,0,0,0,
                         0,0,0,0.02,0.02,
                         1,0,0,0,9.0644e-02,-0.8667e-02,8.2423e-01], dtype=torch.float32, requires_grad=True)
  init_qd = torch.tensor([0,0,0,0,
                          0,0,0,0,0,
                          0,0,0,0,0,0], dtype=torch.float32, requires_grad=False)

  grav = pd.TinyVector3(pd.Utils.zero(), pd.Utils.zero(), pd.Utils.fraction(-1, 980))

  dt = 1/100.
  finer = 1
  world.dt = pd.Utils.scalar_from_double(dt/finer)
  n_step = 1000
  dof_u = 9
  n = dof_u * n_step
  # init_tau = torch.zeros([n], dtype=torch.float32, requires_grad=True)


  dim_control = 7

  param_ctl = torch.normal(mean=0, std=0.03,
    size=(n_step, dim_control), dtype=torch.float32, requires_grad=True)
  # param_ctl = torch.zeros([n_step, dim_control], dtype=torch.float32, requires_grad=True)

  target_joint = [0.2679, 0.2135, 0.6616]
  target_joint = torch.tensor(target_joint, dtype=torch.float32, requires_grad=True)

  world.adj_initialize(grav, n_step, dof_u)
  optimizer = torch.optim.Adam([param_ctl], lr=0.001)
  frameskip_gfx_sync = 16

  # fo = open("ours_throw.txt", "w")
  end_epoch = 60
  for steps in range(500):

    init_tau = construct_ini_tau(param_ctl, n_step)

    sync_counter = frameskip_gfx_sync
    optimizer.zero_grad()
    q, qd = init_q, init_qd
    for sim_step in range(n_step):
      tau = init_tau[sim_step]

      for i in range(finer):
        q, qd = api_diff.sim_layer(q, qd, tau, world)
        clip_dim = 5
        q = torch.cat(
          [q[:clip_dim],
          torch.clamp(q[clip_dim:-7], -0.1, 0.1),
          q[-7:]
          ], axis=0
          )
      sync_counter = sync_counter + 1
      if steps == end_epoch:
        if sync_counter > frameskip_gfx_sync:
          world.sync_visual_meshcat(sim_step)
          sync_counter = 0
    if steps == end_epoch:
      print("recording done")
      exit()
    

    joints = api_diff.get_joints(q, world)
    
    loss = get_loss(joints, target_joint, param_ctl)
    print(steps, " loss =", loss)
    # fo.write( "{}\n".format(loss))

    loss.backward()
    optimizer.step()
  print("done")
  print('tau=',init_tau.reshape([-1,5]))
  print(init_tau.min(), init_tau.max())
  np.save('./data/push/tau.npy', init_tau.reshape([-1,5]).detach().numpy())
  while True:
    continue


if __name__ == "__main__":
  app.run(main)
