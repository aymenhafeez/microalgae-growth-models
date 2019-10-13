import sys
from daetools.pyDAE import *
from time import localtime, strftime
from pyUnits import m, kg, g, s, K, Pa, mol, J, W, kJ, hour, l


class lightModel(daeModel):

    def __init__(self, Name, Parent=None, Description=""):
        daeModel.__init__(self, Name, Parent, Description)
        self.Pin = daeParameter(
            "Pin", g/l, self, "Initial photosynthetic rate")
        self.Pm = daeParameter(
            "Om",   g/l, self, "Maximum photosynthetic rate")
        self.Ik = daeParameter("Ik",   g/l, self, "Light-saturated irradiance")

        self.P = daeVariable("P", molar_concentration_t,
                             self, "Photosynthetic rate")
        self.I = daeVariable("I", unit(), self, "Irradiance")

    def DeclareEquations(self):
        Pin = self.Pin()
        Pm = self.Pm()
        Ik = self.Ik()

        P = self.P()
        I = self.I()

        dP_dI = dI(self.P())

        eq = self.CreateEquation("P", "")
        eq.Residual = dP_dI - Ik * (Pm - P)


class simLightModel(daeSimulation):

    def __init__(self):
        daeSimulation.__init__(self)
        self.m = lightModel("light_model_dae")
        self.m.Description = __doc__

    def SetUpParametersAndDomains(self):
        self.m.Pin.SetValue(0 * g/l)
        self.m.Pm.SetValue(1.654824711 * g/l)
        self.m.Ik.SetValue(34.25 * unit())

    def SetUpVariables(self):
        self.m.P.AssignInitialCondition(0 * g/l)
        self.m.I.AssignInitialCondition(0 * unit())


def run(**kwargs):
    simulation = simLightModel()
    return daeActivity.simulate(simulation, reportingInterval=600,
                                timeHorizon=3*60*60, **kwargs)


if __name__ == "__main__":
    guiRun = False if (
        len(sys.argv) > 1 and sys.argv[1] == 'console') else True
    run(guiRun=guiRun)
