"""calc_peptide.

1.0.0 221104 OKamotoH

"""

import copy
from pathlib import Path
import re
import pandas as pd
from compound import Constants, Unit, Compound


DIR_NAME = Path(__file__).resolve().parent
FILENAME = f'{DIR_NAME}/peptide.xlsx'


class AminoAcid(Unit):
    """アミノ酸の情報のクラス.

    Attributes
    ----------
    letter : str
        アミノ酸の三文字表記
    is_modified_aa : str
         修飾アミノ酸か否かの情報

    See Also
    --------
    Unit: ユニットについて
    Formula : 組成式について

    """

    def __init__(
        self,
        name: str = None,
        note: str = None,
        letter: str = None,
        hydrophilicity: float = 0,
        pKa: float = 0,
        extinction280: int = 0,
        is_hydrophobic: bool = False,
        is_acidic: bool = False,
        is_basic: bool = False,
        is_tfasalt: bool = False,
        **kwargs
    ):
        """初期化メソッド."""
        super().__init__(name=name, note=note, **kwargs)
        self.letter = letter
        self.pKa = pKa
        self.hydrophilicity = hydrophilicity
        self.extinction = extinction280
        self.is_hydrophobic = is_hydrophobic
        self.is_acidic = is_acidic
        self.is_basic = is_basic
        self.is_tfasalt = is_tfasalt
        if len(name) == 1:
            self.is_modified_aa = False
        else:
            self.is_modified_aa = True


class Terminal(Unit):
    """Peptideの末端情報.

    Attributes
    ----------
    **kwargs :
        組成式情報はここで受け取る

    See Also
    --------
    Unit: ユニットについて
    Formula : 組成式について

    """

    def __init__(
        self,
        name: str = None,
        note: str = None,
        **kwargs
    ):
        """初期化メソッド."""
        super().__init__(name=name, note=note, **kwargs)


class ConstantsForPeptide(Constants):
    """ConstatnsをPeptideの分子量計算ができる用に拡張したクラス.

    Attributes
    ----------
    dict_aminoacids : dict[str, AminoAcid]
        アミノ酸一文字表記からアミノ酸のインスタンスを取り出す辞書.
    dict_aminoacids3 : dict[str, AminoAcid]
            アミノ酸一文字表記からアミノ酸のインスタンスを取り出す辞書.
    dict_midifiedaminoacids : dict[str, AminoAcid]
        特殊アミノ酸の名称からアミノ酸のインスタンスを取り出す辞書.
    dict_terminals : dict[str, AminoAcid]
        ペプチドの末端情報のインスタンスを取り出す辞書.

    See Also
    --------
    Constants

    """

    def __init__(self, filename_constants):
        """初期化メソッド."""
        super().__init__(filename_constants)
        sheet_aminoacid = pd.read_excel(
            filename_constants, sheet_name="aminoacid",
            dtype={'is_hydrophobic': bool, 'is_acidic': bool, 'is_basic': bool}
        ).fillna(
            {"extinction280": 0, "note": "",
             "hydrophilicity": 0, "pka": 0,
             "is_hydrophobic": False,
             "is_acidic": False,
             "is_basic": False,
             }
        )
        self.dict_aminoacids_1lett, self.dict_aminoacids_3lett =\
            self.make_dict_aminoacids(sheet_aminoacid)
        self.dict_terminals = self.make_dict_terminals(
            filename_constants)
        sheet_aminoacid.index = sheet_aminoacid["name"]
        self.names_aminoacid = list(sheet_aminoacid.index)
        self.names_acid = list(
            sheet_aminoacid[sheet_aminoacid["is_acidic"]].index)
        self.names_basic = list(
            sheet_aminoacid[sheet_aminoacid["is_basic"]].index)
        self.names_hydrophobic = list(
            sheet_aminoacid[sheet_aminoacid["is_hydrophobic"]].index)

    @staticmethod
    def make_dict_aminoacids(
        sheet_aminoacid,
    ) -> dict[str, AminoAcid]:
        """アミノ酸のインスタンスをまとめた辞書を作成する.

        See Also
        --------
        AminoAcid: アミノ酸のクラス

        """
        sheet_aminoacid["obj"] = sheet_aminoacid.apply(
            lambda x: AminoAcid(**x), axis=1)
        dict_aa_1lett = dict(
            zip(sheet_aminoacid["name"], sheet_aminoacid["obj"]))
        dict_aa_3lett = dict(
            zip(sheet_aminoacid["letter"], sheet_aminoacid["obj"]))
        return dict_aa_1lett, dict_aa_3lett

    @staticmethod
    def make_dict_terminals(filename_constant) -> dict[str, Terminal]:
        """末端情報に関するインスタンスの辞書を作成する.

        See Also
        --------
        Terminal: アミノ酸のクラス

        """
        sheet_aminoacid = pd.read_excel(
            filename_constant, sheet_name="terminal"
        ).fillna("")
        sheet_aminoacid["obj"] = sheet_aminoacid.apply(
            lambda x: AminoAcid(**x), axis=1)
        return dict(zip(sheet_aminoacid["name"], sheet_aminoacid["obj"]))


class SequencePeptide:
    """ペプチドの配列に関するクラス.

    Attributes
    ----------
    constants : ConstantsForPeptide
        ペプチド計算用の定数.
    sequence : str
        ペプチド配列の文字列
    sequence_table : pd.DataFrame
        配列計算用にDataFrame形式に変換したもの.
    cyclo : int
        環を巻いている場合の環の数.
    conv_dict : dict[str, str]
        a,b,c がどの特殊アミノ酸に対応しているか示す辞書.
    length : int
        配列長.

    """

    def __init__(
        self, constants, sequence, sequence_table, cyclo=0
    ):
        """初期化メソッド."""
        self.constants = constants
        self.sequence = sequence
        self.sequence_table = sequence_table
        self.cyclo = cyclo
        self.length = len(self.sequence_table)

    @classmethod
    def create_fm_1letter(cls, sequence, constants, cyclo):
        """ペプチドの文字列の配列（１文字表記）から配列のインスタンスを作成する.

        Parameters
        ----------
        sequence : str
            アミノ酸を一文字で表す形式の配列.
        constants : ConstantsForPeptide
            ペプチド計算用の定数.
        cyclo : int
            環状配列の場合の環の数.
        conv_dict : dict[str, dict]
            特殊ペプチドの変換用の辞書

        """
        sequence = cls.__remove_space(sequence)
        sequence_table = cls.conv_to_table_fm_1letter(constants, sequence)
        return cls(
            constants=constants,
            sequence=sequence,
            sequence_table=sequence_table,
            cyclo=cyclo
        )

    @classmethod
    def create_fm_3letter(cls, sequence, constants, cyclo, sep="-"):
        """ペプチドの文字列の配列（１文字表記）から配列のインスタンスを作成する.

        Parameters
        ----------
        sequence : str
            アミノ酸を一文字で表す形式の配列.
        constants : ConstantsForPeptide
            ペプチド計算用の定数.
        cyclo : int
            環状配列の場合の環の数.
        conv_dict : dict[str, dict]
            特殊ペプチドの変換用の辞書

        """
        sequence = cls.__remove_space(sequence)
        sequence_table = cls.conv_to_table_fm_3letter(
            constants, sequence, sep
        )
        return cls(
            constants=constants,
            sequence=sequence,
            sequence_table=sequence_table,
            cyclo=cyclo
        )

    @staticmethod
    def conv_to_table_fm_1letter(constants, sequence):
        """表形式のSequenceから一文字表記のペプチドに変換する辞書."""
        if re.search(r'\(.+?\)', sequence):
            normal_unit = re.split(r'\(.+?\)', sequence)
            modified_unit = re.findall(r'\((.+?)\)', sequence)
            sequence_list = list()
            for i in range(len(normal_unit)):
                try:
                    sequence_list += list(normal_unit[i])
                    sequence_list += [modified_unit[i]]
                except Exception:
                    if len(normal_unit) < len(modified_unit):
                        sequence_list += list(normal_unit[i])
                    else:
                        pass
        else:
            sequence_list = list(sequence)
        sequece_table = pd.DataFrame(sequence_list, columns=["word"])
        sequece_table["obj"] = sequece_table["word"].apply(
            lambda x: constants.dict_aminoacids_1lett[x]
        )
        return sequece_table[["obj"]]

    @staticmethod
    def conv_to_table_fm_3letter(constants, sequence, sep="-"):
        """表形式のSequenceから一文字表記のペプチドに変換する辞書."""
        if re.search(r'\(.+?\)', sequence):
            normal_unit = re.split(r'\(.+?\)', sequence)
            modified_unit = re.findall(r'\((.+?)\)', sequence)
            sequence_list = list()
            for i in range(len(normal_unit)):
                try:
                    sequence_list += normal_unit[i].split(sep)
                    sequence_list += [modified_unit[i]]
                except Exception:
                    if len(normal_unit) < len(modified_unit):
                        sequence_list += normal_unit[i].split(sep)
                    else:
                        pass
        else:
            sequence_list = sequence.split(sep)
        sequece_table = pd.DataFrame(sequence_list, columns=["word"]).dropna()
        sequece_table = copy.deepcopy(
            sequece_table[sequece_table["word"] != ""])
        sequece_table["obj"] = sequece_table["word"].apply(
            lambda x: constants.dict_aminoacids_3lett[x]
        )
        return sequece_table[["obj"]]

    @staticmethod
    def conv_to_1letter_fm_table(sequence_table):
        """一文字表記のペプチド配列を出力."""
        sequence_table = copy.deepcopy(sequence_table)
        sr_seq = sequence_table["obj"].apply(lambda x: x.name)
        return "".join(sr_seq.to_list())

    @staticmethod
    def conv_to_3letter_fm_table(sequence_table):
        """一文字表記のペプチド配列を出力."""
        sequence_table = copy.deepcopy(sequence_table)
        sr_seq = sequence_table["obj"].apply(lambda x: x.letter)
        return "-".join(sr_seq.to_list())

    def get_formula(self):
        """組成式を得る."""
        return self.sequence_table["obj"].sum().formula

    def get_extinction(self):
        """280nmのモル吸光係数を得る."""
        series_extinction = self.sequence_table["obj"].apply(
            lambda x: x.extinction
        )
        return series_extinction.sum()

    def get_formula(self):
        """組成式を得る."""
        return self.sequence_table["obj"].sum().formula

    def get_extinction(self):
        """280nmのモル吸光係数を得る."""
        series_extinction = self.sequence_table["obj"].apply(
            lambda x: x.extinction
        )
        return series_extinction.sum()

    def get_dict_aminoacid_count(self):
        """アミノ酸の数をカウントした辞書を作成する.

        Return
        ------
        dict[str, int]
            例: {"A": 2, "C": 5, ...}

        """
        aa_word = self.sequence_table["obj"].apply(lambda x: x.name)
        count_aa = pd.Series(
            aa_word.value_counts(),
            index=self.constants.names_aminoacid
        ).fillna(0)
        return count_aa

    def get_tfasalt(self):
        """TFA塩になった際の塩の数を計算する.

        Return
        ------
        int

        """
        aminoacidcount = self.get_dict_aminoacid_count()
        return aminoacidcount[self.constants.names_basic].sum()

    def get_acidic(self):
        """酸性官能基を有するアミノ酸の数を返す.

        Return
        ------
        int

        """
        aminoacidcount = self.get_dict_aminoacid_count()
        return aminoacidcount[self.constants.names_acid].sum()

    def get_basic(self):
        """塩基性官能基を有するアミノ酸の数を返す.

        Return
        ------
        int

        """
        aminoacidcount = self.get_dict_aminoacid_count()
        return aminoacidcount[self.constants.names_basic].sum()

    def get_hydrophobic(self):
        """疎水性アミノ酸の数を計算する.

        Return
        ------
        int

        """
        aminoacidcount = self.get_dict_aminoacid_count()
        return aminoacidcount[self.constants.names_hydrophobic].sum()

    @staticmethod
    def __remove_space(sequence):
        """スペースや改行など余分な物をそぎ落とす."""
        sequence = sequence.replace(" ", "").replace("\n", "")
        sequence = sequence.replace("\t", "")
        return sequence


class Peptide(Compound):
    """ペプチドのクラス.

    Attributes
    ----------
    cyclo : int
        環状ペプチドの場合の環の数
    sequence : SequenceForPeptide
        ペプチド配列のインスタンス
    n_terminal : str
        N末端の情報
    c_terminal : str
        C末端の名称

    See Also
    --------
    compound.Compound : 化合物の基本のクラス.

    """

    def __init__(
        self,
        constants: ConstantsForPeptide,
        sequence: SequencePeptide,
        name: str = "",
        id_sequence: str = "",
        note: str = "",
        deltams: float = 0.0,
        cyclo: int = 0,
        n_terminal: str = "H",
        c_terminal: str = "OH",
        **kwargs
    ):
        """初期化メソッド."""
        self.constants = constants
        self.n_terminal = constants.dict_terminals[n_terminal]
        self.c_terminal = constants.dict_terminals[c_terminal]
        self.cyclo = cyclo
        self.sequence = sequence
        self.extinction = self.sequence.get_extinction()
        formula = self.sequence.get_formula()\
            + self.n_terminal.formula + self.c_terminal.formula
        formula.H -= self.cyclo * 2
        super().__init__(
            formula=formula,
            name=name,
            note=note,
            deltams=deltams,
            constants=self.constants
        )
        self.id_sequence = id_sequence
        self.aminoacids = self.sequence.get_dict_aminoacid_count()
        self.tfasalt = self.calc_tfasalt()
        self.acohsalt = self.calc_acohsalt()
        self.acid = self.sequence.get_acidic()
        self.base = self.sequence.get_basic()
        self.hydrophobic = self.sequence.get_hydrophobic()
        self.ratio_hydrophobic = self.hydrophobic / self.sequence.length
        self.hydrophilicity = self.calc_hydrophilicity()

    def calc_tfasalt(self):
        """TFA塩になる場合の塩の数を返す."""
        salt = self.sequence.get_tfasalt()
        if self.n_terminal.name == "H":
            salt += 1
        return salt

    def calc_acohsalt(self):
        """AcOH塩になる場合の塩の数を返す（仮）."""
        salt = self.aminoacids["K"] + self.aminoacids["R"]
        salt += self.aminoacids["H"]/2
        salt -= (self.aminoacids["D"] + self.aminoacids["E"])
        if self.n_terminal.name == "H":
            salt += 1
        if self.c_terminal.name == "OH":
            salt -= 1
        return salt

    def calc_hydrophilicity(self):
        self.hydrophilicity = 0
        for i in self.aminoacids.index:
            self.hydrophilicity += \
                self.constants.dict_aminoacids_1lett[i].hydrophilicity * \
                self.aminoacids[i]
        self.hydrophilicity = self.hydrophilicity / self.sequence.length
        return self.hydrophilicity

    @classmethod
    def create_fm_1letter(
        cls,
        sequence: str,
        constants: ConstantsForPeptide,
        name: str | None = None,
        id_sequence: str | None = None,
        note: str | None = None,
        deltams: float = 0.0,
        cyclo: int = 0,
        n_terminal: str = "H",
        c_terminal: str = "H",
    ):
        """一文字表記の配列情報からペプチドのインスタンスを作成するメソッド."""
        sequence = SequencePeptide.create_fm_1letter(
            constants=constants,
            sequence=sequence,
            cyclo=cyclo)
        return cls(
            constants=constants,
            name=name,
            note=note,
            deltams=deltams,
            cyclo=cyclo,
            id_sequence=id_sequence,
            sequence=sequence,
            n_terminal=n_terminal,
            c_terminal=c_terminal,
        )

    @classmethod
    def create_fm_3letter(
        cls,
        sequence: str,
        constants: ConstantsForPeptide,
        name: str = None,
        note: str = None,
        deltams: float = 0.0,
        cyclo: int = 0,
        n_terminal: str = "H",
        c_terminal: str = "H",
    ):
        """一文字表記の配列情報からペプチドのインスタンスを作成するメソッド."""
        sequence = SequencePeptide.create_fm_3letter(
            constants=constants,
            sequence=sequence,
            cyclo=cyclo)
        return cls(
            constants=constants,
            name=name,
            note=note,
            deltams=deltams,
            cyclo=cyclo,
            sequence=sequence,
            n_terminal=n_terminal,
            c_terminal=c_terminal,
        )

    def expand_data_for_peptide(self) -> dict:
        """SummarySheet出力用に情報を抜きだす関数."""
        dict_data = {
            "ID": self.id_sequence,
            "Name": self.name,
            "Note": self.note,
            "Length": self.sequence.length,
            "N_terminal": self.n_terminal,
            "Sequence": self.sequence.sequence,
            "C_terminal": self.c_terminal,
            "Cyclo": self.cyclo,
            "Formula": self.formula_str,
            "Ex.MS": self.ex_ms,
            "Mol.W": self.mol_w,
            "ε(280 nm)": self.extinction,
            "TFA": self.tfasalt,
            "AcOH": self.acohsalt,
            "acid": self.acid,
            "base": self.base,
            "hydrophobic": self.hydrophobic,
            "ratio_hydrophobic": self.hydrophobic/self.sequence.length,
            "hydrophilicity": self.hydrophilicity
        }
        dict_data.update(self.sequence.get_dict_aminoacid_count())
        return dict_data

    def expand_data_for_ligand(self) -> dict:
        """SummarySheet出力用に情報を抜きだす関数."""
        dict_data = {
            "Name": self.name,
            "Note": self.note,
            "Formula": self.formula_str,
            "Ex.MS": self.ex_ms,
            "Mol.W": self.mol_w,
        }
        dict_data.update(self.formula.__dict__)
        return dict_data

    def __str__(self):
        """文字列として呼び出された場合の挙動の定義."""
        return f"{self.name} Exact: {self.ex_ms:.4f}; MW:{self.mol_w:.2f}"


if __name__ == '__main__':
    const = ConstantsForPeptide(FILENAME)
    seq1 = SequencePeptide.create_fm_1letter(
        "LGHHIA(K-PEG3-Azide)GKHGY(K-PEG3-Azide)TFG", const, cyclo=0)
    # seq = SequencePeptide.create_fm_3letter("Leu-Ala-(Lys-PEG3-Azide)-Glu", const, cyclo=0)
    pep = Peptide.create_fm_3letter(
        sequence="Leu-Ala-(Lys-PEG3-Azide)-Glu", constants=const)
