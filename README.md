# Reward_Inv_Recursive_Learning
 Use inverse reinforcement learning to determine the value function. 
 Handed in on 2024-05-06 as final project for Cognitive AI. 

# Boltzmann Rationality analysis 
 The script BoltzmannTester accesses all available data records and fits each one to a probability distribution using Bayes' theorem, the beta distribution, and Boltzmann rationality, identifying the beta parameter through several methods. These tasks are accomplished by the getBoltzmann function. Then, it takes the best-fit probability distribution from one recording and tests that on a different recording with the function giveBoltzmann. 

# Apprenticeship IRL in MDP "Decision Space" 
 The script DecisionSpace uses the apprenticeship algorithm to find the value function of binary decisions in one recording session picked at random. 

# Apprenticeship IRL in "Grid World" 
 The script ApprenticeRL uses the apprenticeship algorithm to find the value function of 2D eye movements in one recording session picked at random. 
 Various external helper functions are employed, including: 
  * ReinforcementLearnGrid: implements Q learning provided a value function 
  * getStateSpace: formats recorded eye and cue movement data into the defined continuous state space
  * updateState: update a discrete state with an action
  * gridifyState, ungridState, updateGridState: translate to/from or handle states in the discrete "grid world" space
  * isEqualState: compare states
  * (VisStateSpace: shows states and actions visually)  
