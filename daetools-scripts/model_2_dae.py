import sys
import numpy as np
from daetools.pyDAE import *
from time import localtime, strftime
from pyUnits import m, kg, s, K, Pa, mol, J, W, kJ, hour, l


class model2(daeModel):

    def __init__(self, Name, Parent=None, Description=''):
        daeModel.__init__(self, Name, Parent, Description)
        self.muMax = daeParameter(
            'muMax', day**(-1), self, 'Maximum sepcific growth rate')
        self.ks = daeParameter('ks', mol * m**(-3), self,
                               'Liquid phase mass transfer coefficient')
        self.ki = daeParameter('ki', mol * 1000, self, 'Inhibition constant')
        self.k = daeParameter('k', mol * 1000, self, 'Inhibition constant')
        self.Yi = daeParameter('Yi', mol * kg**(-1), self, 'Yield coefficient')
        self.kla = daeParameter('kla', s**(-1), self,
                                'Mass transfer coefficient')
        self.h = daeParameter('h', Pa * (m**3) * s**(-1),
                              self, 'Henrys constant')
        self.sgin = daeParameter(
            'sgin', Pa, self, 'Initial CO2 partial pressure')
        self.i0 = daeParameter('i0', day**(-1), self, 'Initial irradiance')
        self.a = daeParameter('a', (m**3) * kg**(-1), self, 'Coefficient')

        self.x = daeVariable('x', molar_concentration_t,
                             self, 'Biomass conentration')
        self.sl = daeVariable('sl', molar_concentration_t,
                              self, 'Liquid CO2 conentration')
        self.sg = daeVariable('sg', molar_concentration_t,
                              self, 'CO2 partial pressure')
        self.mu = daeVariable('mu', molar_concentration_t,
                              self, 'Specific growth rate')
        self.i = daeVariable('i', molar_conentration_t, self, 'Irradiance')
        self.t = daeVariable('t', time_1, self, 'Time')

    def DeclareEquations(self):
        muMax = self.muMax()
        ks = self.ks()
        ki = self.ki()
        k = self.k()
        Yi = self.Yi()
        kla = self.kla()
        h = self.h()
        sgin = self.sgin()
        i0 = self.i0()
        a = self.a()

        x = self.x()
        sl = self.sl()
        sg = self.sg()
        mu = self.mu()
        i = self.i()
        t = self.t()

        dx_dt = dt(self.x())
        dsl_dt = dt(self.sl())
        dsg_dt = dt(self.sg())

        eq = self.CreatEquation('mu', '')
        eq.Residual = mu - muMax * sl / \
            ((sl + ks + (sl**2) / ki) * i / (i + k))

        eq = self.CreatEquation('i', '')
        eq.Residual = i - i0 / (a * x * (1 - np.exp(-a * x)))

        eq = self.CreatEquation('x', '')
        eq.Residual = dx_dt - mu * x

        eq = self.CreatEquation('sl', '')
        eq.Residual = dsl_dt - kla * ((sg / h) - sl) - (Yi * dx_dt)

        eq = self.CreatEquation('sg', '')
        eq.Residual = dsg_dt - sgin - kla * ((sg / h) - sl)


class simModel(daeSimulation):

    def __init__(self):
        daeSimulation.__init__(self)
        self.m = model2('model2')
        self.m.Description = __doc__

    def SetUpParametersAndDomains(self):
        self.m.muMax.SetValue(0.5 * 1 / day)
        self.m.ks.SetValue(1 * mol / m**(-3))
        self.m.ki.SetValue(3 * mol * 1000)
        self.m.k.SetValue(14 * mol*1000)
        self.m.Yi.SetValue(0.5 * mol / kg)
        self.m.kla.SetValue(0.00095 * 1 / s)
        self.m.h.SetValue(0.00316 * Pa * m**(3) / s**(-1))
        self.m.sgin.SetValue(0.06 * Pa)
        self.m.i0.SetValue(75 * 1 / day)
        self.m.a.SetValue(0.014 * m**3 / kg)

    def SetUpVariables(self):
        self.m.time_f.AssignValue(1.0 / 16.0)

        self.m.x.SetInitialCondition(0.03 * kg / l)
        self.m.sl.SetInitialCondition(0 * Pa * 1000)
        self.m.sg.SetInitialCondition(17 * Pa * 1000)


def run(**kwargs):
    simulation = simTutorial()
    return daeActivity.simulate(simulation, reportingInterval=600,
                                timeHorizon=3*60*60, **kwargs)


if __name__ == "__main__":
    guiRun = False if (
        len(sys.argv) > 1 and sys.argv[1] == 'console') else True
    run(guiRun=guiRun)
