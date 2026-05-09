

# 装置関係情報
class Column:
  def __init__(self, diameter: int, height: int, fillrate: float = 100):
    self.diameter = diameter
    self.height = height
    self.area = (self.diameter/20) ** 2 * 3.14
    self.volume = self.area * self.height / 10
    self.fillrate = fillrate
    self.bedvolume = self.volume * fillrate / 100
    self.bv = self.bedvolume
    self.cv = self.volume

# 試薬関係情報
class Reagent:
  def __init__(
    self, name: str, mw: float, conc: float, solvent: dict[str, float]
  ):
    self.name = name
    self.mw = mw
    self.conc = conc
    self.conc_unit = "mol/L"
    self.solvent = solvent

  def print_conc(self):
    print(f"{self.conc} mol/L")
    if len(self.solvent) != 1:
      solvs = ":".join(self.solvent.keys())
      ratio = ":".join(map(str, self.solvent.values()))
      solvent = f"({solvs}={ratio})"
    else:
      solvent = f"{list(self.solvent.keys())[0]}"
    print(f" in {solvent}")


class ReagentDeblocking:
  def __init__(self, name: str, conc: float, solvent: dict[str, float]):
    self.name = name
    self.conc = conc
    self.conc_unit = "vol%"
    self.solvent = solvent

  def print_conc(self):
    print(f"{self.conc} vol%")
    if len(self.solvent) != 1:
      solvs = ":".join(self.solvent.keys())
      ratio = ":".join(map(str, self.solvent.values()))
      solvent = f"({solvs}={ratio})"
    else:
      solvent = f"{list(self.solvent.keys())[0]}"
    print(f" in {solvent}")

class Amidite(Reagent):
  def __init__(self, name: str, mw: float, conc: float, solvent: dict[str, float]):
    super().__init__(name, mw, conc, solvent)

class ReagentThio(Reagent):
  def __init__(self, name: str, mw: float, conc: float, solvent: dict[str, float]):
    super().__init__(name, mw, conc, solvent)

class ReagentOxidizer(Reagent):
  def __init__(self, name: str, mw: float, conc: float, solvent: dict[str, float]):
    super().__init__(name, mw, conc, solvent)

class ReagentCapping(Reagent):
  def __init__(self, name: str, mw: float, conc: float, solvent: dict[str, float]):
    super().__init__(name, mw, conc, solvent)


class Constants:
  def __init__(self):
    self.column = Column(diameter=20, height=20, fillrate=80)
    self.amidite = {
      "A": "", "C": "", "G": "", "T": "",
      "a": "", "c": "", "g": "", "u": "",
      "x": "", "y": "", "z": "", "q": "",
    }

#  反応条件
class ParameterDeblocking:
  def __init__(
    self,
    linerflowrate: int, maxflowtime: int, cv_wash: int,
    constants: Constants | None = None
  ):
    if constants == None:
      constants = Constants()
    self.column = constants.column
    self.liner_flowrate = linerflowrate
    self.flowrate = self.flowrate * self.column.area / 60
    self.maxflowtime = maxflowtime
    self.cv_wash = cv_wash
    self.usage = {
      "MeCN": 0.0,
      "Deblocking": 0.0,
    }

  def calc_usage(self):
    self.usage["Deblocking"] += self.flowrate * self.maxflowtime
    self.usage["MeCN"] += self.cv_wash * self.column.cv


class ParameterCoupling:
  def __init__(
    self,
    amidite_reagent: Amidite,
    equiv_amidite: float,
    conc_amidite: float,
    ratio_activator: int,
    time_recycle: float, cv_wash: int,
    constants: Constants | None = None
  ):
    if constants == None:
      constants = Constants()
    self.column = constants.column
    self.flowrate = self.flowrate * self.column.area / 60
    self.cv_wash = cv_wash
    self.usage = {
      "MeCN": 0.0,
      "Deblocking": 0.0,
    }

  def calc_usage(self):
    self.usage["Deblocking"] += self.flowrate * self.maxflowtime
    self.usage["MeCN"] += self.cv_wash * self.column.cv
