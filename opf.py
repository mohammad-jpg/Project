from pyomo.environ import *
import math


#this function prints the results of the solver status of the program
def solverstatus(results):
    if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
        print('feasible')
    elif (results.solver.termination_condition == TerminationCondition.infeasible):
        print('infeasible')
    else:
        print ('Solver Status:',  results.solver.status)

        
        #this function adds 2 lists together and removes the identical elements
def remove_dup(first_list,second_list):
    return(list(dict.fromkeys(first_list + second_list)))   
#returns the limits of the generator
def Pbounds(model,i,t):
     return (0*model.pmin[i]/model.sbase, model.pmax[i]/model.sbase)
#defines the objective function
def objective_function_rule(model):
    return sum(model.b[i]*model.P[i,t] for i in model.i for t in model.t)

def LineBounds(model,lines,bus,node,t):
    
    if (lines,bus,node)  in model.lines:

        return(-model.branchLim[lines,bus,node]/model.sbase,model.branchLim[lines,bus,node]/model.sbase)
    else:
        return Constraint.Skip  
    
def NodesOut_init(model, bus):
    
    retval1 = []
    retval2 = []
    retval = []
    for (l,i,j)  in model.lines:
        if i == bus and j != bus:
            #print(i,j,bus)
            retval1.append(j)
    for (l,j,i)  in model.lines:
        if i == bus and j != bus:
            retval2.append(j)
    retval.append(remove_dup(retval1,retval2))
    #print(str(retval))
    return retval

def Nodeswithgen(model, bus):
    retval = []
    for (gen,node) in model.GBconnect:
        if node == bus:
            retval.append(gen)
    #print(retval) 
    return retval

def linecalc_rule(model,line,bus,node,t):
    if (line,bus,node) in model.lines:
        #print("y")
        return model.flow[line,bus,node,t] == (1/model.branchX[line,bus,node])*(model.delta[bus,t] - model.delta[node,t])
    else:
        return Constraint.Skip
    
def balance_rule(model,bus,t):
    return sum(model.P[i,t] for i in model.GB[bus]) - model.dem[t]*model.Pd[bus]/model.sbase == sum(model.flow[line,b,node,t] for (line,b,node) in model.lines if b == bus)
    

def dcopf(filename):    
        model = AbstractModel()
        model.i = Set(doc='set of Gen units')
        model.b = Param(model.i, doc='b coef of cost')
        model.pmin = Param(model.i, doc='Min gen')
        model.pmax = Param(model.i, doc='Max gen')
        model.t = Set(doc='set of time steps')
        model.dem = Param(model.t, doc='System demand', mutable=True)
        model.GBconnect = Set(dimen=2, doc='Generator connected to bus')

        model.lines = Set(dimen=3,doc = "Lines connecting bus and node")
        model.branchX = Param(model.lines, doc='Max gen')
        model.branchLim = Param(model.lines, doc='Max gen')
        

        model.bus = Set()
        model.node = Set(initialize=model.bus)

        model.sbase = Param(doc='Sbase')
        model.Pd = Param(model.bus, doc='Max gen')
        model.delta = Var(model.bus,model.t, bounds=(-math.pi, math.pi), domain=Reals, doc='angle')


        model.P = Var(model.i,model.t, bounds=Pbounds, domain=Reals)

        model.objective_function = Objective(rule=objective_function_rule, sense = minimize)
        model.GB = Set(model.bus, initialize=Nodeswithgen,doc='Nodes with generator')
  
        model.flow = Var( model.lines,model.t,bounds = LineBounds , domain=Reals, doc='flow')

        model.NodesOut = Set(model.bus, initialize=NodesOut_init,doc='outgoing flow nodes')
        model.constflowcalc = Constraint(model.lines, model.t,rule=linecalc_rule)
        model.constbalance = Constraint(model.bus,model.t , rule=balance_rule)
        opt = SolverFactory('glpk')
        instance = model.create_instance(filename)
        results = opt.solve(instance)
        PowerGen = []
        for i in instance.i:
            for t in instance.t:
                PowerGen.append([i , t, instance.P[i, t].value])


        LineFlow = [] 
        for [l,b,n,t] in instance.flow:
            if instance.flow[l,b,n,t].value >= 0:
                LineFlow.append([l,b,n,t,instance.flow[l,b,n, t].value])

            
        return PowerGen,LineFlow

    
    
