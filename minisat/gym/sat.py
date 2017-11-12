from GymSolver import GymSolver
import numpy as np
import random
from os import listdir
from os.path import isfile, join

class sat(object):
	"""
		this class is a simple wrapper of minisat instance, used in MCTS training as perfect information
	"""
	def __init__(self, sat_dir, max_clause = 100, max_var = 20, mode = 'random'):
		"""
			sat_dir: directory to the sat problems
			max_clause: number of rows for the final state
			max_var: number of columns for the final state
			mode: 'random' => at reset, randomly pick a file from directory
				  'iterate' => at reset, iterate each file one by one
				  'repeat^n' => at reset, give the same problem n times before iterates to the next one
				  'filename' => at reset, repeatly use the given filename
		"""
		print("SAT-v0: at dir {} max_clause {} max_var {} mode {}".format(sat_dir, max_clause, max_var, mode))
		self.sat_dir = sat_dir
		self.sat_files = [join(self.sat_dir, f) for f in listdir(self.sat_dir) if isfile(join(self.sat_dir, f))]
		self.sat_file_num = len(self.sat_files)

		self.max_clause = max_clause
		self.max_var = max_var
		self.observation_space = np.zeros((max_clause, max_var, 2))
		self.action_space = max_var * 2
		self.mode = mode
		if mode.startswith("repeat^"):
			self.repeat_limit = int(mode.split('^')[1])
		elif (mode == "random" or mode == "iterate"): pass
		else:
			try:
				self.file_index = self.sat_files.index(join(self.sat_dir, self.mode))
			except ValueError:
				assert False, "file {} in not in dir {}".format(mode, sat_dir)
		# this class is stateful, by these fields
		self.repeat_counter = 0
		self.iterate_counter = 0

	def reset(self):
		"""
			this function reset the minisat by the rule of mode
		"""
		if self.mode == "random":
			pickfile = self.sat_files[random.randint(0, self.sat_file_num - 1)]
			self.repeat_counter += 1
		elif self.mode == "iterate":
			pickfile = self.sat_files[self.iterate_counter]
			self.iterate_counter += 1
			if self.iterate_counter >= self.sat_file_num:
				self.iterate_counter = 0
				self.repeat_counter += 1
				print("WARNING: iteration of all files in dir {} is done, will restart iteration".format(self.sat_dir))
		elif self.mode.startswith("repeat^"):
			pickfile = self.sat_files[self.iterate_counter]
			self.repeat_counter += 1
			if self.repeat_counter >= self.repeat_limit:
				self.repeat_counter = 0
				self.iterate_counter += 1
				if self.iterate_counter >= self.sat_file_num:
					self.iterate_counter = 0
					print("WARNING: repeated iteration of all files in dir {} is done, will restart iteration".format(self.sat_dir))
		else:
			pickfile = self.sat_files[self.file_index]
			self.repeat_counter += 1
		state = np.zeros((self.max_clause, self.max_var, 2), dtype = np.float32)
		self.S = GymSolver(pickfile)
		self.S.init(np.reshape(state, (self.max_clause*self.max_var*2,)))
		return state
		
		#self.curr_state, self.clause_counter, self.isSolved, self.actionSet = self.parse_state()
		#return state, self.curr_state

	def step(self, decision):
		"""
			this function makes a step based on the parameter input
		"""
		self.S.set_decision(decision)
		state = np.zeros((self.max_clause, self.max_var, 2), dtype = np.float32)
		self.S.step(np.reshape(state, (self.max_clause * self.max_var * 2,)))
		return state, self.S.getReward(), self.S.getDone(), {}
		
	def simulate(self, pi, v):
		"""
			this function makes a simulation step, while providing the pi and v from neural net for the state from the last simulation
			return state (next state to evaluate), bool (need evaluate state, not empty), bool (need more MCTS steps)
		"""
		state = np.zeros((self.max_clause, self.max_var, 2), dtype = np.float32)
		code = self.S.simulate(np.reshape(state, (self.max_clause * self.max_var * 2,)), pi, np.asarray([v], dtype = np.float32))
		if (code == 0): return state, False, False
		if (code == 1): return state, True, False
		if (code == 2): return state, False, True
		if (code == 3): return state, True, True
		assert False, "return code of simulation {} is not one of the designed value".format(code)

	def get_visit_count(self):
		"""
			this function gets the visit count of the root node of MCTS
		"""
		count = np.zeros((self.action_space,), dtype = np.float32)
		self.S.get_visit_count(count)
		return count

"""
def show(array):
	for i in range(100):
		for j in range(20):
			for t in range(2):
				if (array[i, j, t] > 0.5):
					print("{} {} {}".format(i, j, t))
"""

"""
if __name__ == '__main__':
	aa = sat("uf20-91")
	state0 = aa.reset()
	print(state0.shape)
	for i in range(100):
		state, needEnv, needSim = aa.simulate(np.random.randn(40).astype(np.float32), np.float32(0.0))
		print(state.shape)
		print(needEnv)
		print(needSim)
	exit(0)

from sat import sat
from sat import show
import numpy as np 
aa = sat("uf20-91")
state = aa.reset() 
needEnv1 = True 
needSim1 = True
while (needEnv1 or needSim1):
	if not needSim1:
		a, needEnv1, needSim1 = aa.simulate(np.random.randn(40).astype(np.float32), np.float32(0.0))
	else:
		state1, needEnv1, needSim1 = aa.simulate(np.random.randn(40).astype(np.float32), np.float32(0.0))
	print(state1.shape)
	print(needEnv1)
	print(needSim1)


	(aa.step()[0] == state1).sum()
	

	print(state1.shape)
	print(needEnv1)
	print(needSim1)
	state2, needEnv2, needSim2 = aa.simulate(np.random.randn(40).astype(np.float32), np.float32(0.0))
	print(state2.shape)
	print(needEnv2)
	print(needSim2)
	state3, needEnv3, needSim3 = aa.simulate(np.random.randn(40).astype(np.float32), np.float32(0.0))
	print(state3.shape)
	print(needEnv3)
	print(needSim3)
	state4, needEnv4, needSim4 = aa.simulate(np.random.randn(40).astype(np.float32), np.float32(0.0))
	print(state4.shape)
	print(needEnv4)
	print(needSim4)
"""
	


"""
	this function parse the state into sparse matrix (max_clause, max_var, 2) with True for a var in clause 
	NO LONGER NEEDED

def parse_state(self):
	curr_state = np.zeros((self.max_clause, self.max_var, 2), dtype=bool)
	clause_counter = 0
	actionSet = set()
	if not self.S.getDone():
		for line in self.S.getState().split('\n'): # S.getState() gives a string representation of state in cnf format
			if line.startswith("p cnf"):
				header = line.split(" ")
				num_var = int(header[2])
				num_clause = int(header[3])
				assert (num_var <= self.max_var), "num var superseded max var"
			elif line.startswith("c"):
				continue
			elif any(char.isdigit() and (not char == '0') for char in line):
				literals = line.split(" ")
				n = len(literals)
				for j in range(n-1):
					number = int(literals[j])
					nz = 0 if number > 0 else 1
					curr_state[clause_counter, abs(number) - 1, nz] = True
					actionSet.add(number)
				clause_counter += 1
				if clause_counter >= self.max_clause:
					break
	return curr_state, clause_counter, self.S.getDone(), actionSet
"""

"""
if (decision < 0): # this is to say that let minisat pick the decision
	decision = 32767
elif (decision % 2 == 0): # this is to say that pick positive literal
	decision = int(decision / 2 + 1)
else: # this is to say that pick negative literal
	decision = 0 - int(decision / 2 + 1)
if (decision in self.actionSet) or decision == 32767:
	self.S.step(decision)
	self.curr_state, self.clause_counter, self.isSolved, self.actionSet = self.parse_state()
return self.curr_state, self.S.getReward(), self.isSolved, {}
"""