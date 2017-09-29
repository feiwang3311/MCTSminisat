import GymSolver
import numpy as np

class gym_sat(object):
	
	"""
		this class is a gym environment for Reinforcement Learning algorithms
		max_clause: the number of rows in state representation
		max_var: the number of columns in state representation
		satProb: the filename of the SAT problem for this environment
	"""
	def __init__(self, max_clause, max_var, satProb): 
		self.S = GymSolver.GymSolver(satProb) # constructor of GymSolver initialize the SAT problem and out the first state file
		self.satProb = satProb
		self.max_clause = max_clause
		self.max_var = max_var
		self.curr_state, self.clause_counter, self.isSolved, self.actionSet = self.state()

	"""
		this function actually return the current state of the SAT problem
	"""
	def state(self):
		curr_state = np.zeros((self.max_clause, self.max_var), dtype = int)
		clause_counter = 0 # this tracks the current row-to-write (which is permutable!)
		isSolved = False
		actionSet = set() # this set tracks all allowed actions for this state
		with open(self.satProb + "snaps", "r") as file:
			for line in file:
				if line.startswith("p cnf"):
					# this is the header of a cnf problem
					# p cnf 20 90
					header = line.split(" ")
					num_var = int(header[2])
					num_clause = int(header[3])
					if (num_clause == 0 or num_var == 0):
						# sometimes we get degenerated data at the end of the file
						isSolved = True
						break # should break
					assert (num_var <= self.max_var)
					assert (num_clause <= self.max_clause)
				elif line.startswith("c"):
					# other comments line, skip
					continue
				else:
					# clause data line
					# -11 -17 20 0
					literals = line.split(" ")
					n = len(literals)
					for j in range(n-1):
						number = int(literals[j])
						value = 1 if number > 0 else -1
						curr_state[clause_counter, abs(number) - 1] = value
						actionSet.add(number)
					clause_counter += 1
		return curr_state, clause_counter, isSolved, actionSet

	"""
		this function initialize the environment and return the state
	"""
	def initialize(self):
		# TODO: what if self.isSolved is true??
		if (self.isSolved):
			printf("this problem is solved by simplification!\n")
		return self.curr_state, self.isSolved, self.S.getReward(), self.actionSet

	"""
		this function reset the environment and return the initial state
	"""
	def reset(self):
		self.S = GymSolver.GymSolver(self.satProb)
		self.curr_state, self.clause_counter, self.isSolved, self.actionSet = self.state()
		# TODO: what if self.isSolved is true??
		return self.initialize()

	"""
		this function make a step based on parameter input
	"""
	def step(self, decision):
		# TODO: what if self.isSolved is true??
		if decision in self.actionSet:
			self.S.step(decision)
			self.curr_state, self.clause_counter, self.isSolved, self.actionSet = self.state()
			return self.curr_state, self.S.getDone() or self.isSolved, self.S.getReward(), self.actionSet
		# TODO: what if decision is not in actionSet??
		else:
			printf("step decision is not within the action set.\n")	
			return self.initialize() 