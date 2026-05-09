"""PBテーマ用の配列を出力する.

必要ファイル
--------
oligomer.py : オリゴマー配列計算用プログラム
compound.py : 分子量計算用プログラム
constants.xlsx : 糖、塩基、リンカーなどの情報保存用
for_rena.xlsx : 計算管理ファイル

InputFile
---------
for_rena.xlsx : 計算管理ファイル 必須： id_sequence, name, sequence

Output
------
./データベース/data/に oligomer.csv

"""

import datetime
from functools import wraps
import os
from pathlib import Path
import pandas as pd
from oligomer import (
    Base, Sugar, Oligomer, ConstantsForSequence, Sequence, SequenceForOP)


DIR_NAME = Path(__file__).resolve().parent
FILENAME_SEQUENCE = f'{DIR_NAME.parent}/データベース/12_ss_oligonucleotide_forPB.xlsx'
FILENAME_CONSTANT = f'{DIR_NAME}/constants.xlsx'


def stop_watch(func):
    """工程の経過時間出力用Wrapper."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        """時間出力部."""
        print(f"\n>> Start   : {func.__name__:^16} |")
        start = datetime.datetime.now()
        result = func(*args, **kwargs)
        deltatime = datetime.datetime.now() - start
        minitues = int(deltatime.total_seconds()/60)
        sec = deltatime.seconds - 60 * int(deltatime.total_seconds()/60)
        msec = f"{deltatime.microseconds:06d}"[:-4]
        print(f">> End     : {func.__name__:^16} | ",
              f"RAP   : {minitues:02d}′{sec:02d}″{msec:}")
        return result
    return wrapper


class SequenceForRena(Sequence):
    """SequenceのclassをRenaで配列作成時に使用できるように拡張したクラス."""

    @staticmethod
    def check_sugar_and_convert_base(
        sugar_tobe_checked: Sugar,
        base_tobe_checked: Base,
        constants: ConstantsForSequence
    ) -> Base:
        """
        糖の種類から、塩基をTからUへ、CからmCをに変換すべきか判断する関数.
        適宜変換し有るべき塩基のオブジェクトを返す.

        Parameters
        ----------
        sugar_tobe_checked : Sugar
            チェック対象の糖
        base_tobe_checked : Base
            チェック対象の塩基
        constants : ConstantsForSequence
            定数群（塩基オブジェクトの参照用）

        Returns
        -------
        Base
            正しい塩基をとして返す。

        """

        def should_base_be_mecytosine(
            sugar_tobe_checked: Sugar,
            base_tobe_checked: Base,
            constants: ConstantsForSequence
        ) -> bool:
            """MeCytosineとすべきものがCytosineになっている場合 Trueを返す."""
            cytosine = constants.dict_bases["C"]
            should_be_change = False
            if not sugar_tobe_checked.has_cytosine:
                if base_tobe_checked == cytosine:
                    should_be_change = True
            return should_be_change

        def should_base_be_uracil(
            sugar_tobe_checked: Sugar,
            base_tobe_checked: Base,
            constants: ConstantsForSequence
        ) -> bool:
            """UracilとすべきものがThymineになっている場合 Trueを返す."""
            thymine = constants.dict_bases["T"]
            should_be_change = False
            if not sugar_tobe_checked.has_thymine:
                if base_tobe_checked == thymine:
                    should_be_change = True
            return should_be_change

        if should_base_be_mecytosine(
                sugar_tobe_checked, base_tobe_checked, constants
                ):
            checked_base = constants.dict_bases["mC"]
        elif should_base_be_uracil(
                sugar_tobe_checked, base_tobe_checked, constants
                ):
            checked_base = constants.dict_bases["U"]
        else:
            checked_base = base_tobe_checked
        return checked_base

    @classmethod
    def create_fm_baseline_as_gapmer(
        cls,
        sequence_baseline: str,
        sugar_wing: str,
        sugar_gap: str,
        linker_wing: str,
        linker_gap: str,
        wing: list[int],
        constants: ConstantsForSequence,
    ) -> Sequence:
        """配列・修飾の情報からGapmer出力を生成するプログラム.

        Parameters
        ----------
        sequence_baseline : str
            ACGTで表記した塩基文字列
        wing : list[int]
            wingの数
        sugar_wing : str
            wingの糖部
        sugar_gap : str, optional
            gapの糖部 The default is "RNA".
        linker_wing : str, optional
            wing部のリンカー. The default is "PS".
        linker_gap : str, optional
            gap部のリンカー. The default is "PS".
        constants : ConstantsForSequence, optional
            定数群. The default is None.

        Returns
        -------
        obj_sequence : SequenceForRena
            配列情報から作成した配列オブジェクト
        """
        if sum(wing) < len(sequence_baseline):
            # 文字列情報から各種オブジェクトへの変換.
            obj_sugar_wing = constants.dict_sugars[sugar_wing]
            obj_sugar_gap = constants.dict_sugars[sugar_gap]
            obj_link_wing = constants.dict_linkers[linker_wing]
            obj_link_gap = constants.dict_linkers[linker_gap]
            oligomer = cls.create_fm_sequence_ns(sequence_baseline, constants)
            # 糖部の変更
            oligomer.sequence_table["sugar"] = obj_sugar_gap
            oligomer.sequence_table["sugar"][:wing[0]] = obj_sugar_wing
            oligomer.sequence_table["sugar"][-wing[1]:] = obj_sugar_wing
            # Linker部の変更 wingのみ修飾する
            oligomer.sequence_table["linker"][:wing[0]] = obj_link_wing
            oligomer.sequence_table["linker"][
                wing[0]:-wing[1]-1] = obj_link_gap
            oligomer.sequence_table["linker"][-wing[1]-1:-1] = obj_link_wing
            # 塩基部の変更
            for i_mer in oligomer.sequence_table.index:
                i_sugar = oligomer.sequence_table.loc[i_mer, "sugar"]
                i_base = oligomer.sequence_table.loc[i_mer, "base"]
                checked_base = cls.check_sugar_and_convert_base(
                    i_sugar, i_base, constants
                    )
                oligomer.sequence_table.loc[i_mer, "base"] = checked_base
            obj_sequence = cls.create_fm_sequence_table(
                oligomer.sequence_table, constants)
        return obj_sequence

    def is_complementary(
        self,
        complementary_sequence: Sequence,
        print_checksheet: bool = False,
    ) -> bool:
        """対象配列がこの配列と相補的であるかを判断する."""
        # この配列の相補的な配列を 3->5 の表記で作成.
        complementary_baseline = self.get_complementary_baseline(
            rna=True, reverse=True)
        # 確認対象である相補的な配列のbaselineを得る
        complementary_sequence_baseline = complementary_sequence.get_baseline(
            rna=True, reverse=False
            )
        result = complementary_baseline == complementary_sequence_baseline
        if print_checksheet:
            baseline = self.get_baseline(rna=True, reverse=False)
            baseline_comprementary = complementary_sequence.get_baseline(
                rna=True, reverse=True
                )
            print(f"5'-{self.sequence_ns}-3'")
            print(f"5'-{baseline}-3'")
            print(f"   {self.length*'|'}")
            print(f"3'-{baseline_comprementary}-5'")
        return result


class RenaOligomer(Oligomer):
    """管理表for_serviceを処理するためのオリゴのインスタンス.

    Oligomerのクラスを拡張し、IDや名称のAttiributeを追加した。

    Attributes
    ----------
    id_sequence : str
        配列ID
    name_baseline : str
        配列名（リガンド無しの名称）

    See Also
    --------
    Oligomer : オリゴマー

    """

    def __init__(
        self,
        id_sequence: str,
        name_baseline: str,
        constants: ConstantsForSequence,
        sequence: Sequence,
        modification_5: str = "H",
        modification_3: str = "H",
        name: str = None,
        note: str = None,
        deltams: float = 0.0,
        extinction: bool = True,
        sodiumform: bool = False,
        **kwargs,
    ):
        """初期化メソッド."""
        super().__init__(
            constants=constants,
            name=name,
            note=note,
            deltams=deltams,
            sequence=sequence,
            modification_5=modification_5,
            modification_3=modification_3,
            extinction=extinction,
            sodiumform=sodiumform,
            )
        self.id_sequence = id_sequence
        self.name_baseline = name_baseline
        self.sequence_op, self.checksheeet = self.set_opsequence()

    @classmethod
    def read_renafile(
        cls,
        constants: ConstantsForSequence,
        name: str,
        memo: str,
        sequence: str,
        modification_5: str,
        modification_3: str,
        baseline_sequence: str,
        id_sequence: str,
        **kwargs,
    ):
        """Rena管理ファイルからの読み込み."""
        sequence = SequenceForRena.create_fm_sequence_ns(sequence, constants)
        flp = cls(
            name=name,
            note=memo,
            sequence=sequence,
            modification_5=modification_5[:8],
            modification_3=modification_3[:8],
            deltams=0,
            constants=constants,
            extinction=True,
            sodiumform=False,
            name_baseline=baseline_sequence[10:],
            id_sequence=id_sequence,
            )
        flp.id_baseline_sequence = baseline_sequence[:8]
        return flp

    def set_opsequence(self):
        """OPSequenceの作成."""
        sequence = SequenceForOP.create_fm_instance(self.sequence)
        sequence_op = sequence.create_sequence_op_fm_sequence_table()
        checksheet = sequence.make_seq_sequence_opchecksheet()
        return sequence_op, checksheet

    def expand_for_summary(self):
        """サマリーテーブル用に展開."""
        summary_dict = {
            "ID Sequence": self.id_sequence,
            "Name": self.name,
            "memo": self.note,
            "id_mod5": self.modification_5.id_modification,
            "Modification_5": self.modification_5.name,
            "id_baseline": self.id_baseline_sequence,
            "sequence_ns": self.sequence.sequence_ns,
            "id_mod3": self.modification_3.id_modification,
            "Modification_3": self.modification_3.name,
            "Length": self.sequence.length,
            "Formula": self.formula_str,
            "ExactMS": self.ex_ms,
            "Mol.W.": self.mol_w,
            "Extinct": self.extinction,
            "nmol/OD": self.nmol_od,
            "µg/OD": self.ug_od,
            "Sequence(OP)": self.sequence_op,
            "sequence_gd": self.sequence.get_seq_gd(),
            }
        summary_dict.update(self.checksheeet)
        return summary_dict


class BatchSheet:
    """Rena管理シートのクラス."""

    def __init__(self, filename, constants):
        """初期化メソッド."""
        self.constants = constants
        self.sheet = self.read_sheet(filename)
        self.data_for_calc_flp = None
        self.summary_flp = None

    @staticmethod
    def read_sheet(filename):
        """配列管理ファイルの読み込み, Modificationの情報を書き出す."""
        sheet_sequence = pd.read_excel(
            filename, sheet_name="sequence", header=0,
            ).fillna("")
        return sheet_sequence

    @stop_watch
    def calc_flp(self):
        """FLPの計算をしてSumamrySheetを作成する."""
        self.data_for_calc_flp = self.sheet.copy()
        self.data_for_calc_flp["object"] = self.data_for_calc_flp.apply(
            lambda x_flp_info: RenaOligomer.read_renafile(
                constants=self.constants, **x_flp_info),
            axis=1
            )
        self.output_flp_summary()
        return self.sheet

    def output_flp_summary(self):
        """FLPのサマリーを出力する."""
        self.summary_flp = self.data_for_calc_flp.apply(
            lambda x_flp_df: x_flp_df["object"].expand_for_summary(),
            axis=1, result_type="expand")
        dirname = rf"{DIR_NAME.parent}/データベース/data"
        self.summary_flp.to_csv(
            rf"{dirname}/12_ss_oligonucleotide.csv",
            index=False, encoding="utf-8-sig"
            )
        return self.summary_flp

    def calc_all(self):
        """全部計算する用."""
        self.calc_flp()

    def get_obj(self, id_sequence):
        obj = self.data_for_calc_flp[
            self.data_for_calc_flp["id_sequence"] == id_sequence
            ]["object"]
        return list(obj)[0]


def run_all():
    const = ConstantsForSequence(
        filename_constant=FILENAME_CONSTANT,
        filename_sequence=FILENAME_SEQUENCE,
        )
    sheet = BatchSheet(FILENAME_SEQUENCE, const)
    sheet.calc_all()
    # os.startfile(rf"{DIR_NAME.parent}/データベース/オリゴ核酸管理.xlsx")


def get_sodium(seq_id):
    """Na塩の情報を得る."""
    oligo = sheet.get_obj(seq_id)
    sodium_form = oligo.get_sodiumform()
    print(f"{oligo.name=}")
    print(f"{sodium_form.formula=}")
    print(f"{sodium_form.mol_w=:.2f}")
    return str(sodium_form.formula)


if __name__ == '__main__':
    const = ConstantsForSequence(
        filename_constant=FILENAME_CONSTANT,
        filename_sequence=FILENAME_SEQUENCE,
        )
    sheet = BatchSheet(FILENAME_SEQUENCE, const)
    sheet.calc_all()
    os.startfile(rf"{DIR_NAME.parent}/データベース/12_ss_oligonucleotide_forPB.xlsx")
