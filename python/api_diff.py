import torch
import pydiffarti as pd

class SimLayer(torch.autograd.Function):
	@staticmethod
	def forward(ctx, q, qd, tau, world):
		ctx.save_for_backward(q, qd, tau)
		ans = pd.forward_step(q.detach().numpy(), qd.detach().numpy(), tau.detach().numpy(), world)
		ctx.world = world
		q, qd = ans
		q = torch.tensor(q, dtype=torch.float32)
		qd = torch.tensor(qd, dtype=torch.float32)
		ans = q, qd
		return q, qd

	@staticmethod
	def backward(ctx, dldq, dldqd):
		q, qd, tau = ctx.saved_tensors
		ans = pd.backward_step(q.detach().numpy(), qd.detach().numpy(), tau.detach().numpy(), dldq.detach().numpy(), dldqd.detach().numpy(), ctx.world)
		q, qd, tau = ans
		q = torch.tensor(q, dtype=torch.float32)
		qd = torch.tensor(qd, dtype=torch.float32)
		tau = torch.tensor(tau, dtype=torch.float32)
		
		# print('gradtau=',tau)
		return q, qd, tau, None

sim_layer = SimLayer.apply

class GetJoints(torch.autograd.Function):
	@staticmethod
	def forward(ctx, q, world):
		ctx.save_for_backward(q)
		ans = pd.forward_get_joints(q.detach().numpy(), world)
		ctx.world = world
		ans = torch.tensor(ans, dtype=torch.float32)
		return ans

	@staticmethod
	def backward(ctx, dldans):
		q, = ctx.saved_tensors
		ans = pd.backward_get_joints(q.detach().numpy(), dldans.detach().numpy(), ctx.world)
		ans = torch.tensor(ans, dtype=torch.float32)
		return ans, None

get_joints = GetJoints.apply



