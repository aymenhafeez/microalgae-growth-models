import sys
from daetools.pyDAE import *
from time import localtime, strftime
from pyUnits import m, kg, s, K, Pa, mol, J, W, kJ, hour, l


class model1(daemodel):
    def __init__(self, Name, Parent=None, Description=''):
        daeModel.__init__(self, Name, Parent, Description)
        self.Xin = daeParameter('Xin', (kg*1000)/l, self,
                                'Initial Biomass Concentration')
        self.F = daeParameter('F', l/hour, self, 'Flowrate')
        self.V = daeParameter('X1', l, self, 'Culture working volume')
        self.D = daeParameter('D', hour, self, 'Dilution rate')
        self.mu = daeParameter('mu', hour**(-1), self,
                               'Average sepcific growth rate')

        self.X1 = daeVariable('X1', molar_concentration_t,
                              self, 'Biomass Concentration')
        self.t = daeVariable('t', time_1, self, 'Time')

    def DeclareEquations(self):
        Xin = self.Xin()
        F = self.F()
        V = self.V()
        D = self.D()
        mu = self.mu()

        X1 = self.X1()
        t = self.t()

        dX1_dt = dt(self.X1())

        eq = self.CreateEquation('X1', '')
        eq.Residual = dX1_dt - X * (mu - (F / V)) + (F / V) * Xin


class simModel(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = model1('model_1')
        self.m.Description = __doc__

    def SetUpParametersAndDomains(self):
        self.m.Xin.SetValue(0.028745092 * (kg*1000)/l)
        self.m.F.SetValue(0.00045 * l/hour)
        self.m.V.SetValue(0.15 * l)
        self.m.mu.SetValue(0.01321570192 * 1/hour)

    def SetUpVariables(self):
        self.m.time_f.AssignValue(1.0 / 16.0)

        self.m.X.SetInitialCondition(0.028745092 * (kg*1000)/l)


def run(**kwargs):
    simulation = simTutorial()
    return daeActivity.simulate(simulation, reportingInterval=600,
                                timeHorizon=3*60*60, **kwargs)
    availableSimulations = {}
    availableSimulations['model_1'] = loaderFunction
    daeSimulationWebService.runSimulationsAsWebService(availableSimulations)


if __name__ == "__main__":
    guiRun = False if (
        len(sys.argv) > 1 and sys.argv[1] == 'console') else True
    run(guiRun=guiRun)
