"""組成式計算用の組成式や、化合物、化合物の構成要素のクラス.

20220509 OK-A-moto,H.

"""
import datetime
from pathlib import Path
import pandas as pd

# pylint: disable=C0103

DIR_NAME = Path(__file__).resolve().parent
FILENAME_CONSTANT = f'{DIR_NAME}/constants.xlsx'


class Constants:
  """定数類をまとめたもの.

  最低限の情報を持つ定数群, Compoundなど配列情報を持たないMS計算をするときに使用.

  Attributes
  ----------
  starttime : datetime.datetime
      開始時間
  proton_mol_w : float
      分子量計算用のprotonの質量
  proton_exact : float
      精密質量計算用のprotonの質量
  weight_mol_w : Weight
      分子量計算用の原子の質量をまとめたインスタンス
  weight_exact : Weight
      精密質量計算用の原子の質量をまとめたインスタンス

  """

  def __init__(
      self,
      filename_constant: str = FILENAME_CONSTANT,
  ):
    """初期化メソッド.

    Parameters
    ----------
    filename_constant : str optional
        定数をまとめたExcelファイル名
        Default: FILENAME_CONSTANT
    filename_sequence : str optional
        配列、Ligandなど使用者によって異なる情報をまとめたExcelファイル名
        Default: FILENAME_SEQUENCE

    """
    # time
    self.starttime = datetime.datetime.now()
    # Weight
    self.proton_mol_w = 1.0079
    self.proton_exact = 1.00727645232093
    self.weight_mol_w = Weight.set_value(
        filename_constant=filename_constant,
        mol_w=True
    )
    self.weight_exact = Weight.set_value(
        filename_constant=filename_constant,
        mol_w=False
    )


class Weight:
  """質量の情報のクラス."""

  def __init__(
      self,
      C: float = 0,
      H: float = 0,
      B: float = 0,
      N: float = 0,
      O: float = 0,
      F: float = 0,
      Na: float = 0,
      Si: float = 0,
      P: float = 0,
      S: float = 0,
      Cl: float = 0,
      K: float = 0,
      **kwargs
  ):
    self.C = C
    self.H = H
    self.B = B
    self.Cl = Cl
    self.F = F
    self.K = K
    self.N = N
    self.Na = Na
    self.O = O
    self.P = P
    self.S = S
    self.Si = Si

  @classmethod
  def set_value(
      cls,
      mol_w: bool,
      filename_constant: str = FILENAME_CONSTANT,
  ):
    """値の設定方法を選択しインスタンスを作成する.

    Parameters
    ----------
    mol_w : bool, optional
        精密質量か分子量かを選択する。
        True: 分子量, False: 精密質量
        The default is True.
    filename_constant : str, optional
        値に指定があるときは指定する.
        無い場合は標準の値を使用する
        The default is FILENAME_CONSTANT.
    TYPE
        DESCRIPTION.

    Returns
    -------
    Weight

    """
    sheet_weight = pd.read_excel(
        filename_constant, sheet_name="weight", header=0, index_col='name'
    ).fillna("")
    if mol_w:
      weight = cls(**sheet_weight.loc["mol_w"])
    else:
      weight = cls(**sheet_weight.loc["exact"])
    return weight


class Formula:
  """組成式のクラス."""

  def __init__(
      self,
      C: int = 0,
      H: int = 0,
      B: int = 0,
      Cl: int = 0,
      F: int = 0,
      K: int = 0,
      N: int = 0,
      Na: int = 0,
      O: int = 0,
      P: int = 0,
      S: int = 0,
      Si: int = 0,
      **kwargs
  ):
    self.C = int(C)
    self.H = int(H)
    self.B = int(B)
    self.Cl = int(Cl)
    self.F = int(F)
    self.K = int(K)
    self.N = int(N)
    self.Na = int(Na)
    self.O = int(O)
    self.P = int(P)
    self.S = int(S)
    self.Si = int(Si)

  def get_ms(
      self,
      mol_w=True,
      constants=None
  ) -> float:
    """組成式から求められる質量を計算するメソッド.

    Parameters
    ----------
    mol_w : bool
      計算する質量を分子量か精密質量か選択する。
      The default is True
      True: 分子量(Mol.Weightを計算する)
      False: 精密質量(Exact MSを計算する)
    constants: Constants
      定数を束ねたクラスConstants ここからWeightの基準値を取得するために指定
      The default is None
      Noneの時は標準の質量（defaultのfile_constantの値）を使用

    Return
    ------
    calcdms : float
        Formulaの原子の個数に対してそれぞれ原子のWeightをかけた和

    """
    if constants is None:
      weight = Weight.set_value(mol_w=mol_w)
    else:
      if mol_w is True:
        weight = constants.weight_mol_w
      else:
        weight = constants.weight_exact
    calcdms = self * weight
    return calcdms

  def __str__(self):
    """文字列として呼び出された場合組成式の文字列として値を返す."""
    formula_str = ''
    # 原子と数値を書き出していく。
    for i_atom, i_number in self.__dict__.items():
      atom = i_atom
      number = i_number
      if number == 0:
        formula_str += ""
      elif number == 1:
        formula_str += f"{atom}"
      else:
        formula_str += f'{atom}{int(number)}'
    return formula_str

  def __add__(self, other):
    """組成式通しの足し算を行う。同じ原子同士で足し合わせる."""
    if isinstance(other, self.__class__):
      C = self.C + other.C
      H = self.H + other.H
      B = self.B + other.B
      Cl = self.Cl + other.Cl
      F = self.F + other.F
      K = self.K + other.K
      N = self.N + other.N
      Na = self.Na + other.Na
      O = self.O + other.O
      P = self.P + other.P
      S = self.S + other.S
      Si = self.Si + other.Si
      formula = Formula(C, H, B, Cl, F, K, N, Na, O, P, S, Si)
    elif hasattr(other, "formula"):
      formula = self + other.formula
    else:
      raise NotImplementedError()
    return formula

  def __sub__(self, other):
    """組成式通しの引き算を行う。同じ原子同士で足し合わせる."""
    if isinstance(other, self.__class__):
      C = self.C - other.C
      H = self.H - other.H
      B = self.B - other.B
      Cl = self.Cl - other.Cl
      F = self.F - other.F
      K = self.K - other.K
      N = self.N - other.N
      Na = self.Na - other.Na
      O = self.O - other.O
      P = self.P - other.P
      S = self.S - other.S
      Si = self.Si - other.Si
      formula = Formula(C, H, B, Cl, F, K, N, Na, O, P, S, Si)
    else:
      raise NotImplementedError()
    return formula

  def __mul__(self, other):
    """組成式のかけ算.

    >> かけ算の相手がWeightの場合
    質量計算に使用.各原子に対して各原子のWeightを書けて質量を返す。
    >> かけ算の相手がintの場合
    各原子に対してint倍した組成式を返す。

    """
    if isinstance(other, Weight):
      weight = 0
      weight += self.C * other.C
      weight += self.H * other.H
      weight += self.B * other.B
      weight += self.Cl * other.Cl
      weight += self.F * other.F
      weight += self.K * other.K
      weight += self.N * other.N
      weight += self.Na * other.Na
      weight += self.O * other.O
      weight += self.P * other.P
      weight += self.S * other.S
      weight += self.Si * other.Si
      return_value = weight
    elif isinstance(other, int):
      C = self.C * other
      H = self.H * other
      B = self.B * other
      Cl = self.Cl * other
      F = self.F * other
      K = self.K * other
      N = self.N * other
      Na = self.Na * other
      O = self.O * other
      P = self.P * other
      S = self.S * other
      Si = self.Si * other
      return_value = Formula(C, H, B, Cl, F, K, N, Na, O, P, S, Si)
    else:
      return_value = None
      raise NotImplementedError()
    return return_value

  def __rmul__(self, other):
    """組成式のかけ算.

    >> かけ算の相手がWeightの場合
    質量計算に使用.各原子に対して各原子のWeightを書けて質量を返す。
    >> かけ算の相手がintの場合
    各原子に対してint倍した組成式を返す。

    """
    if isinstance(other, Weight):
      weight = 0
      weight += self.C * other.C
      weight += self.H * other.H
      weight += self.B * other.B
      weight += self.Cl * other.Cl
      weight += self.F * other.F
      weight += self.K * other.K
      weight += self.N * other.N
      weight += self.Na * other.Na
      weight += self.O * other.O
      weight += self.P * other.P
      weight += self.S * other.S
      weight += self.Si * other.Si
      return_value = weight
    elif isinstance(other, int):
      C = self.C * other
      H = self.H * other
      B = self.B * other
      Cl = self.Cl * other
      F = self.F * other
      K = self.K * other
      N = self.N * other
      Na = self.Na * other
      O = self.O * other
      P = self.P * other
      S = self.S * other
      Si = self.Si * other
      return_value = Formula(C, H, B, Cl, F, K, N, Na, O, P, S, Si)
    else:
      return_value = None
      raise NotImplementedError()
    return return_value

  def __neg__(self):
    """符号反転用."""
    C = -self.C
    H = -self.H
    B = -self.B
    Cl = -self.Cl
    F = -self.F
    K = -self.K
    N = -self.N
    Na = -self.Na
    O = -self.O
    P = -self.P
    S = -self.S
    Si = -self.Si
    return Formula(C, H, B, Cl, F, K, N, Na, O, P, S, Si)


class Unit:
  """物質を構成する要素のクラス.

  Attributes
  ----------
  name : str
      ユニットの名前
  formula : Formula
      ユニットの組成式のインスタンス
  note : str
      ユニットの補足情報

  See Also
  --------
  Formula: 組成式について

  """

  def __init__(
      self,
      name: str,
      note: str,
      **kwargs
  ):
    """初期化メソッド."""
    self.name = name
    self.note = note
    self.formula = Formula(**kwargs)

  def __str__(self):
    """文字列として呼ばれたとき."""
    return self.name

  def __add__(self, other):
    """下記情報を持つUnitとして返す.

    Note
    ----
    主にDeletionを作成するとき、欠損するヌクレオチドを足していく
    name: other.name, self.name
    note: other.note, self.note
    formula: self.formula + other.formula

    """
    if isinstance(other, Unit):
      name = f"{self.name}, {other.name}"
      note = f"{self.note}, {other.note}"
      formula = self.formula + other.formula
      return Unit(name, note, **formula.__dict__)
    else:
      raise NotImplementedError()


class Compound:
  """化合物のクラス.

  Attributes
  ----------
  formula : Formula
      化合物の組成式のインスタンス
  formula_str : str
      組成式の文字列
  name : str
      化合物の名称
  note : str
      化合物の説明
  deltams : float
      組成不明の質量差
  constants : Constants
      定数のインスタンス

  Note
  ----
  組成式を持つもの。
  OligomerはCompoundに配列情報などを拡張したもの。

  """

  def __init__(
      self,
      formula: Formula,
      name: str = None,
      note: str = None,
      deltams: float = 0.0,
      constants=None,
  ):
    """初期化メソッド.

    Parameters
    ----------
    formula : Formula
        化合物の組成式.
    name : str, optional
        化合物の名前. The default is None.
    note : str, optional
        化合物の説明. The default is None.
    deltams : float, optional
        組成不明の質量差. 標準は0.0
    constants : Constants, optional
        定数を束ねた物、質量の基準値抽出用に設定.
        標準はNoneでそのときにはConstantsのインスタンスを作成する。

    """
    self.constants = constants
    self.formula = formula
    self.formula_str = self.get_formula_str(formula, deltams)
    self.name = name
    self.note = note
    self.deltams = deltams
    self.mol_w = formula.get_ms(
        mol_w=True, constants=self.constants
    )
    self.ex_ms = formula.get_ms(
        mol_w=False, constants=self.constants
    )

  @staticmethod
  def get_formula_str(formula, deltams=0):
    """deltamsを加味したformula_strを出力する."""
    if deltams > 0:
      deltams_str = f"+{deltams:.4f}"
    elif deltams < 0:
      deltams_str = f"-{deltams:.4f}"
    else:
      deltams_str = ""
    return f"{formula}{deltams_str}"

  def make_mz_dict(
      self,
      num_z: int = 1,
      mol_w: bool = True,
  ) -> dict():
    """m/zの値を展開する.

    Parameters
    ----------
    z_num : int
        何価までm/zの値を算出するか
    mol_w : bool
        m/zの値を分子量or精密質量で算出するかの判断
        True: 分子量, False: 精密質量
    constants : Constants
        算出する際の定数のインスタンス
        The default is None:
        Noneのときメソッド中で標準のインスタンスを作成する.

    Return
    ------
    mz_list : dict[int, float]
        z価までm/zの値を計算したList

    """
    if mol_w is True:
      ms_value = self.mol_w
      proton = self.constants.proton_mol_w
    else:
      ms_value = self.ex_ms
      proton = self.constants.proton_exact
    num_z = max(num_z, 1)  # 1 merの時計算されないので対応.
    mz_table = pd.DataFrame(range(-1, -(num_z + 1), -1), columns=["z"])
    mz_table["m/z"] = mz_table["z"].apply(
        lambda x: (ms_value + x * proton) / -x
    )
    return dict(zip(mz_table["z"], mz_table["m/z"]))

  def expand_data_for_summary(
      self,
      num_z=1,
      mz_data: bool = False,
      mol_w: bool = True,
  ) -> dict:
    """SummarySheet出力用に情報を抜きだす関数."""
    dict_data = {
        "Name": self.name,
        "Note": self.note,
        "Formula": self.formula_str,
        "ExactMS": self.ex_ms,
        "Mol.W": self.mol_w,
    }
    if mz_data is True:
      dict_mz = self.make_mz_dict(
          num_z=num_z,
          mol_w=mol_w,
      )
      dict_data.update(dict_mz)
    return dict_data


if __name__ == '__main__':
  # ForDEBUG
  test_const = Constants()
  test_formula = Formula(C=2, H=6, O=1)
  aa = test_formula.get_ms(mol_w=True)
  bb = Compound(formula=test_formula, constants=test_const)
  cc = bb.make_mz_dict(num_z=5)
  ttt = Formula(C=2, H=4, O=1)
  aaaa = Formula(C=160, H=320, O=80)
