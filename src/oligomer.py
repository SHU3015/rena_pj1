"""Oligomerの基本のクラス.

Base, Sugar, Linker, Modification, Nucleotideのクラス.
Sequence, SequenceForSort, SequenceForOP, Oigomerのクラス.
FLPの計算までする事が可能.

20220509 OK-A-moto,H.
20230314 OligoPilot100 CheckSheet変更

"""
from __future__ import annotations
import copy
import re
import sys
from pathlib import Path
import pandas as pd
from compound import Constants, Formula, Unit, Compound


DIR_NAME = Path(__file__).resolve().parent
FILENAME_SEQUENCE = f"{DIR_NAME.parent}/データベース/02_ss_oligonucleotide.xlsx"
FILENAME_CONSTANT = Path(f"{DIR_NAME}/constants.xlsx").resolve()


class Base(Unit):
    """塩基の情報のクラス.

    Attribute
    ---------
    extinct : str
        塩基のモル吸光係数

    See Also
    --------
    Unit: ユニットについて
    Formula : 組成式について

    """

    def __init__(
        self, name: str | None = None, note: str | None = None, extinct: int = 0, **kwargs
    ):
        """初期化メソッド."""
        super().__init__(name=name, note=note, **kwargs)
        self.extinct = extinct


class Sugar(Unit):
    """糖の情報のクラス.

    Attributes
    ----------
    has_bridge : bool
        架橋型核酸か否か. MSMSの付加物判断に使用
    sign_ns : str
        糖名のns配列における１文字での略称
    singn_gd : str
        糖名のgd配列における１文字での略称
    **kwargs :
        組成式情報はここで受け取る

    See Also
    --------
    Unit: ユニットについて
    Formula : 組成式について

    """

    def __init__(
        self,
        name: str | None = None,
        note: str | None = None,
        is_bridge: bool | None = False,
        sign_ns: str | None = None,
        sign_gd: str | None = None,
        is_cytosine: bool | None = None,
        is_thymine: bool | None = None,
        **kwargs,
    ):
        """初期化メソッド.

        Parameters
        ----------
        name : str
            糖の名称
        note : str
            糖の補足
        extinct : int
            塩基のモル吸光係数
        has_bridge : bool
            架橋型核酸か否か. MSMSの付加物判断に使用
        sign_ns : str
            糖名のns配列における１文字での略称
        singn_gd : str
            糖名のgd配列における１文字での略称

        """
        super().__init__(name=name, note=note, **kwargs)
        self.sign_ns = sign_ns
        self.sign_gd = sign_gd
        self.has_cytosine = is_cytosine
        self.has_thymine = is_thymine
        if is_bridge == 1:
            self.has_bridge = True
        else:
            self.has_bridge = False


class Linker(Unit):
    """Linkerの情報のクラス.

    Attributes
    ----------
    mark : str
        リンカーの略称記号
        例 PS: ^, PB: * など
    **kwargs :
        組成式情報はここで受け取る

    See Also
    --------
    Unit: ユニットについて
    Formula : 組成式について

    """

    def __init__(self, name: str | None = None, note: str | None = None, mark: str = "", **kwargs):
        """初期化メソッド."""
        super().__init__(name=name, note=note, **kwargs)
        self.mark = mark


class Modification(Unit):
    """塩基の情報のクラス.

    Attribute
    ---------
    location : str
        5'の付加物か3'の付加物かを記す。
    type : str
        msms付加物か、Ligandかを記す。
    salt : int
        塩になり得る活性プロトンが何個あるかを示す。
    **kwargs :
        組成式情報はここで受け取る

    See Also
    --------
    Unit: ユニットについて
    Formula : 組成式について

    """

    def __init__(
        self,
        name: str | None = None,
        note: str | None = None,
        id_modification: str | None = None,
        location: int | None = None,
        type_modification: str | None = None,
        salt: int = 0,
        extinction: int = 0,
        note_for_ns: str | None = None,
        note_for_gd: str | None = None,
        **kwargs,
    ):
        """初期化メソッド."""
        super().__init__(name=name, note=note, **kwargs)
        self.id_modification = id_modification
        self.location = location
        self.type = type_modification
        self.salt = salt
        self.extinction = extinction
        self.note_for_ns = note_for_ns
        self.note_for_gd = note_for_gd


class Nucleotide(Unit):
    """ヌクレオチド情報.

    Attributes
    ----------
    name : str
        A, G(L)^などNS-Rena表記
    note : str
        DNA-A(PO), LNA-G(PS)などの表記
    formula : Formula
        組成式のインスタンス
    base : Base
        構成要素 Baseのインスタンス
    sugar : Sugar
        構成要素 Sugarのインスタンス
    linker : Linker
        構成要素 Linkerのインスタンス

    """

    def __init__(
        self,
        base: Base,
        sugar: Sugar,
        linker: Linker,
    ):
        """初期化メソッド."""
        super().__init__(
            name=self.make_nucleotide_name(base, sugar, linker),
            note=f"{sugar.name}-{base.name}({linker.name})",
        )
        self.formula = base.formula + sugar.formula + linker.formula
        self.base = base
        self.sugar = sugar
        self.linker = linker
        self.name_ns = self.name
        self.name_gd = self.make_nucleotide_name(base, sugar, linker, "gd")

    @staticmethod
    def make_nucleotide_name(base, sugar, linker, name_cmo="ns"):
        """sequence_tableの列（Base, Sugar, Linker）をNucleotideに変換する."""
        if name_cmo == "ns":
            if sugar.name == "RNA":
                basename = base.name.lower()
            else:
                basename = base.name
        elif name_cmo == "gd":
            if base.name == "mC":
                basename = "5"
            elif sugar.name == "DNA":
                basename = base.name.lower()
            else:
                basename = base.name
        sign = sugar.__dict__[f"sign_{name_cmo}"]
        if sign == "":
            nucleotide = f"{basename}{linker.mark}"
        else:
            nucleotide = rf"{basename}({sign}){linker.mark}"
        return nucleotide

    def __str__(self):
        """文字列として呼ばれたとき."""
        return self.name

    def __add__(self, other):
        """ヌクレオチド同士の足し算."""
        if isinstance(other, Unit):
            # 主にDeletionを作成するときの仕様. Unitとして返す.
            name = f"{self.name}, {other.name}"
            note = f"{self.note}, {other.note}"
            formula = self.formula + other.formula
            return Unit(name, note, **formula.__dict__)
        elif isinstance(other, Formula):
            # Formulaとの計算は計算後Formulaのみを返す.
            return self.formula + other
        else:
            print(type(other))
            raise NotImplementedError()


class ConstantsForSequence(Constants):
    """定数群をオリゴ配列の分子量計算をできる用に拡張した.

    Attributes
    ----------
    dict_bases : dict[str, Base]
        塩基の情報のインスタンスが格納された辞書
    dict_sugars : dict[str, Sugar]
        糖の情報のインスタンスが格納された辞書
    dict_linkers : dict[str, Linker]
        Linkerの情報のインスタンスが格納された辞書
    dict_modifications5 : dict[str, Modification]
        5'末に付加しうるリガンドやMSMS付加物の情報のインスタンスが格納された辞書
    dict_modifications3 : dict[str, Modification]
        3'末に付加しうるリガンドやMSMS付加物の情報のインスタンスが格納された辞書
    dict_sign_to_sugar : dict[str, dict[str, str]]
        会社によって異なる糖の一文字の略号から糖名に変換する辞書をまとめたもの
    dict_impurities : dict[ste, pd.DataFrame]
        不純物ファイル名から不純物（組成式差）表を呼び出す辞書

    SeeAlso
    -------
    Constants

    """

    def __init__(
        self,
        filename_constant: str = FILENAME_CONSTANT,
        filename_sequence: str = FILENAME_SEQUENCE,
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
        super().__init__(filename_constant)
        # Unit Dictionarys
        self.dict_bases = self.make_dict_bases(filename_constant)
        self.dict_sugars = self.make_dict_sugars(filename_constant)
        self.dict_linkers = self.make_dict_linkers(filename_constant)
        self.dict_modifications5 = self.make_dict_modification(
            filename_sequence, filename_constant, 5
        )
        self.dict_modifications3 = self.make_dict_modification(
            filename_sequence, filename_constant, 3
        )
        # Conv_dict
        self.dict_sign_to_sugar = dict()
        self.dict_sign_to_sugar["ns"] = self.make_conv_dict(self.dict_sugars, name_cmo="ns")
        self.dict_sign_to_sugar["gd"] = self.make_conv_dict(self.dict_sugars, name_cmo="gd")

    @staticmethod
    def make_dict_bases(filename_constant) -> dict[str, Base]:
        """塩基のインスタンスを作成し辞書を作成する.

        See Also
        --------
        Base: 塩基のクラス

        """
        bases = dict()
        sheet_base = pd.read_excel(filename_constant, sheet_name="Base", header=0).fillna("")
        for i_line in sheet_base.index:
            name = sheet_base.loc[i_line, "name"]
            bases.update({name: Base(**sheet_base.loc[i_line])})
        return bases

    @staticmethod
    def make_dict_sugars(filename_constant) -> dict[str, Sugar]:
        """糖のインスタンスを作成し辞書を作成する.

        See Also
        --------
        Sugar: 糖のクラス

        """
        sugars = dict()
        sheet_sugar = pd.read_excel(filename_constant, sheet_name="Sugar", header=0).fillna("")
        for i_line in sheet_sugar.index:
            name = sheet_sugar.loc[i_line, "name"]
            sugars.update({name: Sugar(**sheet_sugar.loc[i_line])})
        return sugars

    @staticmethod
    def make_dict_linkers(filename_constant) -> dict[str, Linker]:
        """Linkerのインスタンスを作成し辞書を作成する.

        See Also
        --------
        Linker: Linkerのクラス

        """
        linkers = dict()
        sheet_linker = pd.read_excel(filename_constant, sheet_name="Linker", header=0).fillna("")
        for i_line in sheet_linker.index:
            name = sheet_linker.loc[i_line, "name"]
            linkers.update({name: Linker(**sheet_linker.loc[i_line])})
        return linkers

    @staticmethod
    def make_dict_modification(
        filename_sequence, filename_constant, place
    ) -> dict[str, Modification]:
        """リガンドやMSMS付加物の情報のインスタンスが格納された辞書を作成する.

        Parameters
        ----------
        filename_constant : str
            定数をまとめたExcelファイル名
        filename_sequence : str
            配列、Ligandなど使用者によって異なる情報をまとめたExcelファイル名
        place : str
            5'末か3'末かを指定.

        Returns
        -------
        dict_modifications : dict[str, Modification]
            リガンドやMSMS付加物の情報のインスタンスが格納された辞書

        """
        dbname = {"5": "DB04_Modification5", "3": "DB05_Modification3"}
        mod = pd.read_excel(filename_sequence, sheet_name=f'{dbname[f"{place}"]}', header=0).fillna(
            0
        )
        dict_modifications = dict()
        for i_mod in mod.index:
            modification_unit = Modification(
                location=place, type_modification="modification", **mod.loc[i_mod]
            )
            dict_modifications.update({modification_unit.id_modification: modification_unit})
        msms = pd.read_excel(filename_constant, sheet_name=f"msms{place:.0f}", header=0)
        for i_msms in msms.index:
            modification_unit = Modification(
                location=place, type_modification="msms", **msms.loc[i_msms]
            )
            dict_modifications.update({modification_unit.name: modification_unit})
        return dict_modifications

    @staticmethod
    def make_conv_dict(sugars, name_cmo="ns") -> dict[str, str]:
        """糖の一文字の略号から糖名に変換する辞書を作成する.

        Parameters
        ----------
        sugars : dict[str, Sugar]
            dict_sugars
        name_cmo : str, optional
            どのメーカーの辞書を作成するか。 The default is "ns".

        Returns
        -------
        dict_sign_to_sugar : dict[str, str]
            糖の一文字の略号から糖名に変換する辞書

        """
        sugseries = pd.Series(sugars)
        sugseries2 = sugseries.apply(lambda x: x.__dict__[f"sign_{name_cmo}"]).dropna()
        dict_sign_to_sugar = dict(zip(sugseries2, sugseries2.index))
        return dict_sign_to_sugar


class Sequence:
    """塩基配列のクラス.

    Attributes
    ----------
    sequence_ns : str
        日本触媒形式の表記で記された配列情報
    sequence_table : pd.DataFrame
        | 塩基配列を5'から順に１行ずつ 塩基、糖、Linkerの情報に切り分けた表.
        | * index: 5'末端から数えて何塩基目か
        | * column: [base, sugar, linker]
    sequence_nucleotide : pd.DataFrame
        | 塩基配列を5'から順に１行ずつ ヌクレオチドの情報をならべたもの.
        | * index: 5'末端から数えて何塩基目か
        | * columns: ["nucleotide", "name_nucleotide", "note_nucleotide"]
    length : int
        配列長

    """

    def __init__(
        self,
        constants: ConstantsForSequence,
        sequence_ns: str,
        sequence_table: pd.DataFrame,
        sequence_nucleotide: pd.DataFrame,
    ):
        """初期化メソッド.

        Parameters
        ----------
        constants : ConstantsForSequence or ConstantsForCalculation
            定数情報のインスタンス.
        sequence_ns : str
            日本触媒形式の表記で記された配列情報
        sequence_table : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ 塩基、糖、Linkerの情報に切り分けた表.
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["base", "sugar", "linker"]
        sequence_nucleotide : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ ヌクレオチドの情報をならべたもの
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["nucleotide", "name_nucleotide", "note_nucleotide"]

        """
        self.constants = constants
        self.sequence_ns = sequence_ns.strip()
        self.sequence_table = sequence_table
        self.sequence_nucleotide = sequence_nucleotide
        self.length = len(self.sequence_table)

    @classmethod
    def create_fm_sequence_ns(
        cls,
        sequence_ns: str,
        constants: ConstantsForSequence | None = None,
    ):
        """SequenceのインスタンスをNS形式の配列(sequence_ns)から生成する.

        Parameters
        ----------
        sequence_ns : str
            NS形式での配列情報.
        constants : ConstantsForSequence or ConstantsForCalculation
            定数群のインスタンス

        See Also
        --------
        Sequence.conv_to_table_fm_str : 文字列からsequence_tableに変換する
        Sequence.conv_to_nucleotide_fm_table :
            sequence_tableからsequence_nucleotideに変換する

        """
        if constants is None:
            constants = ConstantsForSequence()
        else:
            pass
        sequence_ns = sequence_ns.strip()
        sequence_table = cls.conv_to_table_fm_str(
            constants=constants, sequence=sequence_ns, name_cmo="ns"
        )
        sequence_nucleotide = cls.conv_to_nucleotide_fm_table(sequence_table)
        return cls(
            constants=constants,
            sequence_ns=sequence_ns,
            sequence_table=sequence_table,
            sequence_nucleotide=sequence_nucleotide,
        )

    @classmethod
    def create_fm_sequence_gd(
        cls,
        sequence_gd: str,
        constants: ConstantsForSequence | None = None,
    ):
        """SequenceのインスタンスをGD形式の配列(sequence_gd)から生成する.

        Parameters
        ----------
        sequence_gd : str
            NS形式での配列情報.
        constants : ConstantsForSequence or ConstantsForCalculation
            定数群のインスタンス

        See Also
        --------
        Sequence.conv_to_table_fm_str : 文字列からsequence_tableに変換する
        Sequence.conv_to_nucleotide_fm_table :
            sequence_tableからsequence_nucleotideに変換する

        """
        if constants is None:
            constants = ConstantsForSequence()
        else:
            pass
        sequence_gd = sequence_gd.strip()
        sequence_table = cls.conv_to_table_fm_str(
            constants=constants, sequence=sequence_gd, name_cmo="gd"
        )
        sequence_nucleotide = cls.conv_to_nucleotide_fm_table(sequence_table)
        sequence_ns = cls.conv_to_str_fm_nucleotide(sequence_nucleotide)
        return cls(
            constants=constants,
            sequence_ns=sequence_ns,
            sequence_table=sequence_table,
            sequence_nucleotide=sequence_nucleotide,
        )

    @classmethod
    def create_fm_sequence_table(
        cls,
        sequence_table: pd.DataFrame,
        constants: ConstantsForSequence,
    ):
        """sequence_tableから各種配列情報を作成しインスタンスを生成させる.

        Parameters
        ----------
        sequence_table : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ 塩基、糖、Linkerの情報に切り分けた表.
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["base", "sugar", "linker"]
        constants : ConstantsForSequence or ConstantsForCalculation
            定数群のインスタンス

        See Also
        --------
        Sequence.conv_to_str_fm_nucleotide : sequence_tableから文字列に変換する
        Sequence.conv_to_nucleotide_fm_table :
            sequence_tableからsequence_nucleotideに変換する

        """
        sequence_nucleotide = cls.conv_to_nucleotide_fm_table(sequence_table)
        sequence_ns = cls.conv_to_str_fm_nucleotide(sequence_nucleotide)
        return cls(
            constants=constants,
            sequence_ns=sequence_ns,
            sequence_table=sequence_table,
            sequence_nucleotide=sequence_nucleotide,
        )

    @classmethod
    def create_fm_sequence_nucleotide(
        cls,
        sequence_nucleotide: pd.DataFrame,
        constants: ConstantsForSequence,
    ):
        """sequence_nucleotideから各種配列情報を作成しインスタンスを生成させる.

        Parameters
        ----------
        sequence_nucleotide : pd.DataFrame
            塩基配列を5'から順に１行ずつ ヌクレオチドの情報をならべたもの
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["nucleotide", "name_nucleotide", "note_nucleotide"]
        constants : ConstantsForSequence or ConstantsForCalculation
            定数群のインスタンス

        See Also
        --------
        Sequence.conv_to_str_fm_nucleotide : sequence_tableから文字列に変換する
        Sequence.conv_to_nucleotide_fm_table :
            sequence_tableからsequence_nucleotideに変換する

        """
        sequence_table = cls.conv_to_table_fm_nucleotide(sequence_nucleotide)
        sequence_ns = cls.conv_to_str_fm_nucleotide(sequence_nucleotide)
        return cls(
            constants=constants,
            sequence_ns=sequence_ns,
            sequence_table=sequence_table,
            sequence_nucleotide=sequence_nucleotide,
        )

    @staticmethod
    def conv_to_table_fm_str(
        constants: ConstantsForSequence,
        sequence: str,
        name_cmo: str = "ns",
    ) -> pd.DataFrame:
        r"""sequence_tableに配列（文字列）から変換するメソッド.

        Parameters
        ----------
        sequence : str
            配列情報.
        name_cmo : str
            どのcmo表記の配列かcmo名を選択する。現状はNS、GDに対応。
        constants : Constants
            定数群のインスタンス、conv_dictの取得用

        Returns
        -------
        sequence_table : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ 塩基、糖、Linkerの情報に切り分けた表.
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["base", "sugar", "linker"]

        Note
        ----
        * 塩基の判断
            mCの判別のために、m, 5かどうかのみのチェック。大文字小文字の判別はしていない。
        * 糖の判別
            糖の判別について
                | 塩基の文字列の次の文字列が"()"で有るかどうかを判断する。
                | ⇒（）がある：括弧内の文字を修飾核酸として認識する。
                | ⇒（）がない：塩基が大文字ならDNA、小文字ならRNAとする。
            "()"の認識について、
                塩基の文字のパターンは mC or Xなので、r'[a-zA-Z]{1,2}\(\D\)'として
                正規表現で抽出しようとしていたが、cG(L)の様な配列もヒットしてしまうので、
                mCかどうかを判別した後に、r'[a-zA-Z0-9]\(\D\)'で抽出するようにしている。
        """

        def check_nucleotide(sequence, num):
            """配列文字列左端のヌクレオチドの情報を抜きだす.

            IndexError時は3末のオリゴと判断しLinkerをNoneとする。

            Parameters
            ----------
            sequence : str
                元配列
            num : int
                元配列の何文字目からチェックするかの数

            Returns
            -------
            basename : str
                塩基の名前 check_baseで判断
            sugar_name : str
                糖の名前 check_sugarで判断
            linker_name : str
                リンカーの名前 check_linkerで判断
            skip : int
                塩基、糖、リンカーの情報に使用した文字数。

            """

            def check_base(
                sequence: str,
                num: int,
            ) -> tuple[str, int]:
                """塩基の情報を判断するメソッド.

                mCかどうかを判断。m, 5で無ければ文字をそのまま大文字にして塩基と判断する。

                Returns
                -------
                base_name : str
                    塩基の名前。例 A, mCなど
                skip : int
                    糖の情報に何文字使っているかの値。mC: 2, その他 1

                Note
                ----
                塩基を2文字にする場合、大文字と小文字で異なる塩基を意味する様なパターンは不適.

                """
                if sequence[num] == "m":
                    base_name = "mC"
                    skip = 2
                elif sequence[num] == "5":
                    base_name = "mC"
                    skip = 1
                else:
                    base_name = sequence[num].upper()
                    skip = 1
                return base_name, skip

            def check_sugar(
                sequence: str,
                num: int,
                skip: int,
            ) -> tuple[str, int]:
                r"""糖の情報の判断するメソッド.

                Parameters
                ----------
                num : int
                    配列のどこから判断するか。
                skip : int
                    塩基情報が何文字分であるかの値. check_baseで判断した値.

                Returns
                -------
                sugar_name : str
                    糖の名称
                skip : int
                    塩基と糖の情報に使用した文字数

                Note
                ----
                アルゴリズムについて、
                    | 塩基の文字列の次の文字列が"()"で有るかどうかを判断する。
                    | ⇒（）がある：括弧内の文字を修飾核酸として認識する。
                    | ⇒（）がない：塩基が大文字ならDNA、小文字ならRNAとする。
                "()"の認識について、
                    塩基の文字のパターンは mC or Xなので、r'[a-zA-Z]{1,2}\(\D\)'として
                    正規表現で抽出しようとしていたが、cG(L)の様な配列もヒットしてしまうので、
                    mCかどうかを判別した後に、r'\((.+?)\)'で抽出するようにしている。
                    22/10/17 修飾核酸の複数文字対応した。

                """
                if sequence[num] == "m":
                    checkpoint = num + 1
                else:
                    checkpoint = num
                try:
                    if re.match(r"[a-zA-Z0-9]\(.+?\)", sequence[checkpoint:]):
                        sugar_signs = re.findall(r"\((.+?)\)", sequence[checkpoint:])
                        sugar_sign = sugar_signs[0]
                        sugar_name = constants.dict_sign_to_sugar[name_cmo][sugar_sign]
                        skip += len(sugar_sign) + 2
                    else:
                        if name_cmo == "ns":
                            if sequence[checkpoint].islower() is True:
                                sugar_name = "RNA"
                            elif sequence[checkpoint].isupper() is True:
                                sugar_name = "DNA"
                            else:
                                raise NameError("Sugar名が対象外です。")
                        elif name_cmo == "gd":
                            if sequence[checkpoint].islower() is True:
                                sugar_name = "DNA"
                            elif sequence[checkpoint].isupper() is True:
                                sugar_name = "RNA"
                            else:
                                print(sequence[checkpoint])
                                sys.exit("Error! def CheckSugar")
                        skip += 0
                except IndexError:  # 不要?
                    print(r"この表示が出たら岡本まで連絡ください。@def check_sugar")
                    if sequence[checkpoint].islower() is True:
                        sugar_name = "RNA"
                    elif sequence[checkpoint].isupper() is True:
                        sugar_name = "DNA"
                        skip += 0
                    else:
                        sys.exit("Error! def CheckSugar")
                return sugar_name, skip

            def check_linker(sequence, num, skip) -> tuple[str, int]:
                """Linker情報を判断する.

                Parameters
                ----------
                sequence : str
                    配列の文字列
                num : int
                    何文字目から判断するかの値
                skip : int
                    塩基と糖の情報に使用した文字数

                Returns
                -------
                linker_name : str
                    リンカーの名前
                skip : int
                    塩基、糖、リンカーの情報に使用した文字数

                """
                next_str = sequence[num + skip]
                skip += 1
                if next_str == r"^":
                    linker_name = "PS"
                elif next_str == r"*":
                    linker_name = "PB"
                elif next_str == r"~":
                    linker_name = "PNMe2"
                else:
                    linker_name = "PO"
                    skip -= 1
                return linker_name, skip

            try:
                basename, skip = check_base(sequence, num)
                sugar_name, skip = check_sugar(sequence, num, skip)
                linker_name, skip = check_linker(sequence, num, skip)
            except IndexError:
                basename, skip = check_base(sequence, num)
                sugar_name, skip = check_sugar(sequence, num, skip)
                linker_name = "none"
            return basename, sugar_name, linker_name, skip

        # 5'末端から順に何塩基目のヌクレオチドがどの塩基、糖、リンカーから構成されるか切り分ける。
        list_sequence = list()
        append_list_sequence = list_sequence.append
        num = 0
        while num < len(sequence):
            base, sugar, linker, skip = check_nucleotide(sequence, num)
            append_list_sequence(
                {
                    "base": constants.dict_bases[base],
                    "sugar": constants.dict_sugars[sugar],
                    "linker": constants.dict_linkers[linker],
                }
            )
            num += skip
        sequence_table = pd.DataFrame(list_sequence)
        return sequence_table

    @staticmethod
    def conv_to_str_fm_nucleotide(sequence_nucleotide: pd.DataFrame, name_cmo: str = "ns") -> str:
        """各CMO表記の配列に表形式の配列から変換するメソッド.

        Parameters
        ----------
        sequence_nucleotide : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ ヌクレオチドの情報をならべたもの
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["nucleotide", "name_nucleotide", "note_nucleotide"]
        name_cmo : str, optional
            どのCMOの配列か. The default is "ns".

        Returns
        -------
        sequence_str : str
            各CMO表記の配列の文字列

        """
        sr_nucleotide_str = sequence_nucleotide["nucleotide"].apply(
            lambda x_nucleotide: x_nucleotide.__dict__[f"name_{name_cmo}"]
        )
        return "".join(sr_nucleotide_str.to_list())

    @staticmethod
    def conv_to_nucleotide_fm_table(sequence_table: pd.DataFrame) -> pd.DataFrame:
        """sequence_nucleotideにsequence_tableから変換するメソッド.

        Parameters
        ----------
        sequence_table : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ 塩基、糖、Linkerの情報に切り分けた表.
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["base", "sugar", "linker"]

        Return
        ------
        sequence_nucleotide : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ ヌクレオチドの情報をならべたもの
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["nucleotide", "name_nucleotide", "note_nucleotide"]

        """
        sequence_nucleotide = pd.DataFrame(
            columns=["nucleotide", "name_nucleotide", "note_nucleotide"]
        )
        sequence_nucleotide["nucleotide"] = sequence_table.apply(
            lambda x: Nucleotide(x["base"], x["sugar"], x["linker"]), axis=1
        )
        sequence_nucleotide["name_nucleotide"] = sequence_nucleotide["nucleotide"].apply(
            lambda x: x.name
        )
        sequence_nucleotide["note_nucleotide"] = sequence_nucleotide["nucleotide"].apply(
            lambda x: x.note
        )
        return sequence_nucleotide

    @staticmethod
    def conv_to_table_fm_nucleotide(sequence_nucleotide: pd.DataFrame) -> pd.DataFrame:
        """sequence_tableにsequence_nucleotideから変換するメソッド.

        Parameters
        ----------
        sequence_nucleotide : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ ヌクレオチドの情報をならべたもの
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["nucleotide", "name_nucleotide", "note_nucleotide"]

        Return
        ------
        sequence_table : pd.DataFrame
            | 塩基配列を5'から順に１行ずつ 塩基、糖、Linkerの情報に切り分けた表.
            | * index: 5'末端から数えて何塩基目か
            | * columns: ["base", "sugar", "linker"]

        """

        def expand_to_table(x_line) -> pd.Series:
            """Nucleotideのインスタンスが持つSugar, Base, Linkerの情報を展開するメソッド.

            Parameters
            ----------
            x_line : pd.Series
                5末から数えてn番目のヌクレオチドのインスタンス、名前、略称のシリーズ.

            Returns
            -------
            dict_to_table : pd.Series
                5末から数えてn番目のヌクレオチドのインスタンスをSugar, Base, Linkerのインスタンスに
                展開したSeries.

            """
            nucleotide = x_line["nucleotide"]
            dict_to_table = {
                "base": nucleotide.base,
                "sugar": nucleotide.sugar,
                "linker": nucleotide.linker,
            }
            return dict_to_table

        sequence_table = sequence_nucleotide.apply(expand_to_table, axis=1, result_type="expand")
        return sequence_table

    @staticmethod
    def conv_to_ns_fm_gd(sequence_gd: str, constants: ConstantsForSequence = None) -> str:
        """NS形式の配列表記へ、GD形式の配列表記から変換を行う.

        Parameters
        ----------
        sequence_gd : str
            GD形式の配列文字列

        Returns
        -------
        sequence_ns : str
            NS形式の配列文字列

        See Also
        --------
        Sequence.create_fm_sequence_gd :
            gd形式の配列文字列からインスタンスを作成するクラスメソッド.
            sequence_ns

        """
        if constants is None:
            constants = ConstantsForSequence()
        seq = Sequence.create_fm_sequence_gd(sequence_gd, constants)
        return seq.sequence_ns

    @staticmethod
    def conv_to_gd_fm_ns(sequence_ns: str, constants: ConstantsForSequence = None) -> str:
        """GD形式の配列表記へ、NS形式の配列表記から変換を行う.

        Parameters
        ----------
        sequence_ns : str
            NS形式の配列文字列

        Returns
        -------
        sequence_gd : str
            GD形式の配列文字列


        See Also
        --------
        Sequence.create_fm_sequence_ns :
            gd形式の配列文字列からインスタンスを作成するクラスメソッド.
        Sequence.get_seq_gd :
            GD表記の配列情報を取得する

        """
        if constants is None:
            constants = ConstantsForSequence()
        seq = Sequence.create_fm_sequence_ns(sequence_ns, constants)
        return seq.get_seq_gd()

    def get_seq_reverse(self, is_seq_ns=True):
        """3->5表記の配列を得る.

        Parameters
        ----------
        is_seq_ns : bool
            True NS表記の配列を返す False GD表記の配列を返す

        Returns:
            is_seq_ns TrueならNS表記の配列を返す
            is_seq_ns FalseならGD表記の配列を返す

        """
        seq = copy.deepcopy(self.sequence_table)
        seq["linker"][0] = self.sequence_table["linker"][len(seq) - 1]
        for i in seq.index:
            seq["linker"][i + 1] = self.sequence_table["linker"][i]
        rev_seq = seq[::-1]
        revseq = Sequence.create_fm_sequence_table(rev_seq, self.constants)
        if is_seq_ns:
            return revseq.sequence_ns
        else:
            return revseq.get_seq_gd()

    def get_seq_gd(self) -> str:
        """GD表記の配列情報を取得するメソッド.

        sequence_gdはAttributesに含まれないの取得する関数を用意した。

        See Also
        --------
        Sequence.conv_to_str_fm_nucleotide :
            Sequenceインスタンスのsequence_nucleotideから文字列に変換する関数.

        """
        return Sequence.conv_to_str_fm_nucleotide(self.sequence_nucleotide, "gd")

    def get_nucleotidelist(self) -> list[str]:
        """５末から順に格納されたヌクレオチド名のリストを得るメソッド."""
        return self.sequence_nucleotide["name_nucleotide"]

    def get_formula(self) -> Formula:
        """配列（末端情報無し）の組成式を得るメソッド.

        Returns
        -------
        self.sequence_nucleotide.sum().formula : Formula
            配列（末端情報無し）の組成式

        """
        return self.sequence_nucleotide["nucleotide"].sum().formula

    def get_extinction(self) -> int:
        """配列のモル吸光係数を得るメソッド.

        Returns
        -------
        extinction : int
            この配列のモル吸光係数

        Note
        ----
        各塩基のモル吸光係数の総和に0.9をかけた値.

        """
        sr_extinction = self.sequence_table["base"].apply(lambda x: x.extinct)
        extinction = sr_extinction.sum() * 0.9
        return extinction

    def get_dict_nucleotidecount(self) -> dict[str, int]:
        """ヌクレオチドの個数をカウントした辞書を得るメソッド.

        Return
        ------
        dict[str, int]
            例: {"A^": 2, "C": 5, ...}

        Note
        ----
        sort_seq_nucleotide_and_table後に実行すると整列されてよろしい.

        """
        nucleotide_count = self.sequence_nucleotide["name_nucleotide"].value_counts(sort=False)
        return nucleotide_count

    def get_dict_basecount(self) -> dict[str, int]:
        """塩基数(糖骨格を無視)をカウントした辞書を得るメソッド.

        Return
        ------
        base_count : dict[str, int]
            例: {"A": 2, "C": 5, ...}

        """
        base_count = self.sequence_table["base"].value_counts(sort=True)
        return base_count

    def get_str_basecount(self) -> dict[str, int]:
        """塩基数(糖骨格を無視)をカウントした辞書を得るメソッド.

        Return
        ------
        base_count : dict[str, int]
            例: {"A": 2, "C": 5, ...}

        """
        base_count = self.sequence_table["base"].value_counts(sort=True)
        print(base_count.to_dict())
        str_basecount = ""
        for i in self.constants.dict_bases.keys():
            if i != "none":
                str_basecount += f"{i}: {base_count.get(i)}; "
        return str_basecount

    def get_baseline(self, rna=False, reverse=False) -> str:
        """塩基配列（All-DNA表記）を得るメソッド.

        Parameter
        ---------
        rna : bool
            ACGUで返す場合はTrue, ACGTで返す場合はFalse, defaultはFalse
        reverse : bool
            5->3表記 False, 3->5表記 True

        Returns
        -------
        Baseline : str
            塩基配列

        """

        def conversion(base):
            if rna:
                conv_dict = {"mC": "C", "T": "U"}
            else:
                conv_dict = {"mC": "C", "U": "T"}
            if str(base) in conv_dict.keys():
                base = conv_dict[str(base)]
            return str(base)

        baseline_series = self.sequence_table["base"].apply(conversion)
        baseline = "".join(baseline_series.to_list())
        if reverse:
            baseline = baseline[::-1]
        else:
            baseline = baseline
        return baseline

    def get_complementary_baseline(self, rna=False, reverse=False) -> str:
        """相補鎖のACGTのみの配列表記(Baseline)を得る.

        Parameters
        ----------
        rna : bool
            ACGUで返す場合はTrue, ACGTで返す場合はFalse, defaultはFalse
        reverse : bool
            5->3表記 False, 3->5表記 True

        Returns
        -------
        str
        相補鎖のBaseline

        """
        if rna:
            dict_complementary = {"A": "U", "C": "G", "G": "C", "T": "A", "U": "A"}
        else:
            dict_complementary = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}
        baseline = self.get_baseline(reverse=False)
        complementary = baseline.translate(str.maketrans(dict_complementary))
        return complementary[::-1]


class SequenceForSort(Sequence):
    """配列のソート用のクラス."""

    def __init__(self, constants, sequence_ns, sequence_table, sequence_nucleotide):
        """初期化メソッド."""
        super().__init__(constants, sequence_ns, sequence_table, sequence_nucleotide)
        self.saved_sequence_ns = sequence_ns

    @classmethod
    def create_fm_instance(cls, sequence):
        """SequenceのインスタンスをSequenceForSortに変換する."""
        seq = copy.deepcopy(sequence)
        return cls(seq.constants, seq.sequence_ns, seq.sequence_table, seq.sequence_nucleotide)

    def set_first_coupling_linkage(self, linker=None) -> pd.DataFrame:
        """不純物やDeletion作成用に初回カップリング時のLiker情報を加えた表に変更する.

        Note
        ----
        UnyLinkerに対する初回カップリングの酸化・硫黄化の条件。
        指定が無い場合はNoneとし、自動で判別する。
        配列内にPO結合があれば初回カップリングはPOと見なしそれ意外は"PS"。
        基本Noneで処理。指定する場合はプログラムで変更するのでご連絡ください。

        """
        if linker is None:
            linkername = self.sequence_table["linker"].apply(lambda x: x.name)
            if "PO" in linkername.to_list():
                self.sequence_table["linker"][-1:] = self.constants.dict_linkers["PO"]
            else:
                self.sequence_table["linker"][-1:] = self.constants.dict_linkers["PS"]
        else:
            self.sequence_table["linker"][-1:] = self.constants.dict_linkers[linker]
        self.sequence_nucleotide = self.conv_to_nucleotide_fm_table(self.sequence_table)
        return self

    def sort_sequence(self, linker=None):
        """Deletion作成用にヌクレオチドの並びをソートしたSequenceのインスタンスを作成する."""
        self.set_first_coupling_linkage(linker)
        self.sequence_nucleotide = self.sequence_nucleotide.sort_values("note_nucleotide")
        self.sequence_nucleotide = self.sequence_nucleotide.reset_index(drop=True)
        self.sequence_table = self.conv_to_table_fm_nucleotide(self.sequence_nucleotide)
        return self

    def get_unique_nucleotide(self, linker=None) -> pd.DataFrame:
        """オリゴ中に含まれるヌクレオチドの周囲の列挙.

        Return
        ------
        pd.DataFrame
            ColumnName : values
            object : ヌクレオチドのインスタンス
            note_nucleotide : A^ 表記のヌクレオチド名
            name_nucleotide : DNA-A(PS)表記のヌクレオチド名

        """
        self.set_first_coupling_linkage(linker)
        self.sort_sequence(linker=linker)
        df_unique_nucleotide = copy.deepcopy(self.sequence_nucleotide)
        df_unique_nucleotide = df_unique_nucleotide.drop_duplicates(subset="name_nucleotide")
        return df_unique_nucleotide

    def get_dict_nucleotides(self, linker=None) -> dict[str, Nucleotide]:
        """ヌクレオチドの名称からヌクレオチドのインスタンスを返す辞書.

        Return
        ------
        dict[str, Nucleotide]
            key: A^表記のヌクレオチド名, value: ヌクレオチドのインスタンス

        """
        df_unique_nucleotide = self.get_unique_nucleotide(linker)
        dict_nucleotides = dict(
            zip(df_unique_nucleotide["name_nucleotide"], df_unique_nucleotide["nucleotide"])
        )
        return dict_nucleotides

    def reset_sequence(self):
        """もと配列を復元する."""
        return Sequence.create_fm_sequence_ns(self.constants, self.saved_sequence_ns)


class SequenceForOP(Sequence):
    """SequenceのクラスをOP配列変換用できる様に拡張したクラス."""

    @classmethod
    def create_fm_instance(cls, sequence):
        """SequenceのインスタンスをSequenceForSortに変換する."""
        return cls(
            sequence.constants,
            sequence.sequence_ns,
            sequence.sequence_table.copy(),
            sequence.sequence_nucleotide.copy(),
        )

    def create_sequence_op_fm_sequence_table(self, linker=None) -> str:
        """OP合成用の配列名を作成するプログラム.

        Parameters
        ----------
        sequence_table : pd.DataFrame
            テーブル形式の配列情報

        Returns
        -------
        sequence_op : str
            OPFormatでの配列情報

        Note
        ----
        ・糖の種類が３種類の時は何も返さない。
        ・糖の種類が２種類の時はDNA or RNAをACGT, 修飾核酸をxyzqとする。
        ・糖の種類がDNAとRNAの時はDNAをACGT, RNAをxyzqとする.
        ・POが含まれる場合は末端PO.PSのみならPSとする。

        """

        def convert(line, major_sugar):
            conv_dict_major = {
                "A": "A",
                "C": "C",
                "G": "G",
                "T": "T",
                "U": "T",
                "mC": "C",
            }
            conv_dict_minor = {
                "A": "x",
                "C": "y",
                "G": "z",
                "T": "q",
                "U": "q",
                "mC": "y",
            }
            if line["sugar"] in major_sugar:
                base = conv_dict_major[line["base"].name]
            else:
                base = conv_dict_minor[line["base"].name]
            return f"{base}{line['linker'].mark}"

        seq_for_sort = SequenceForSort.create_fm_instance(self)
        seq_for_sort.set_first_coupling_linkage(linker)
        seq_table = seq_for_sort.sequence_table
        seq_table["sugar"] = seq_table["sugar"].apply(str)
        if len(seq_table["sugar"].unique()) <= 2:
            # 糖の種類を確認し変換
            if ("RNA" and "DNA") in seq_table["sugar"].unique():
                major_sugar = "DNA"
            else:
                major_sugar = ("RNA", "DNA")
            sequence_series = seq_table.apply(lambda x: convert(x, major_sugar), axis=1)
            sequence_op = "".join(sequence_series.to_list())
        else:
            sequence_op = ""
        return sequence_op

    def make_seq_sequence_opchecksheet(self, linker=None) -> dict[str, str]:
        """OPでの合成に使うBlockNameのシリーズを返す.

        Returns
        -------
        sequence_series : dict[str, str]
            OPでの合成に使うBlockNameのシリーズを返す

        """

        def convert(line, major_sugar):
            conv_dict_major = {
                "A": "A",
                "C": "C",
                "G": "G",
                "T": "T",
                "U": "T",
                "mC": "C",
            }
            conv_dict_minor = {
                "A": "x",
                "C": "y",
                "G": "z",
                "T": "q",
                "U": "q",
                "mC": "y",
            }
            # Base
            if line["sugar"] in major_sugar:
                base = conv_dict_major[line["base"].name]
            else:
                base = conv_dict_minor[line["base"].name]
            # Sugar
            if line["sugar"] == "DNA":
                sugar = "DNA"
            elif line["sugar"] == "RNA":
                sugar = "RNA"
                if base == "T":
                    base = "U"
            else:
                sugar = "RNA"
            # Linker
            if str(line["linker"]) == "PS":
                linker = "s"
            else:
                linker = ""
            # data
            return f"{linker}{sugar}-{base}"

        seq_for_sort = SequenceForSort.create_fm_instance(self)
        seq_for_sort.set_first_coupling_linkage(linker)
        seq_table = seq_for_sort.sequence_table[::-1]
        seq_table = seq_table.reset_index(drop=True)
        seq_table["sugar"] = seq_table["sugar"].apply(str)
        if len(seq_table["sugar"].unique()) <= 2:
            # 糖の種類を確認し変換
            if ("RNA" and "DNA") in seq_table["sugar"].unique():
                major_sugar = "DNA"
            else:
                major_sugar = ("RNA", "DNA")
            sequence_series = seq_table.apply(lambda x: convert(x, major_sugar), axis=1)
            sequence_series[0] += "_UNY"
            df_checkblock = pd.DataFrame(sequence_series, columns=["block"])
            df_checkblock["num"] = df_checkblock.index.map(lambda x: f"Block-{x+1:02}")
            dict_checksheet = dict(zip(df_checkblock["num"], df_checkblock["block"]))
        else:
            dict_checksheet = dict()
        return dict_checksheet


class Oligomer(Compound):
    """オリゴマーのクラス.

    At-tributes
    ---------
    constants : constants
        定数群
    name : str
        オリゴの名前
    note : str
        オリゴの補足情報
    sequence : Sequence
        オリゴの配列情報をまとめたインスタンス
    modification_5 : Modification
        5'末端に付加する修飾の情報のインスタンス
    modification_3 : Modification
        3'末端に付加する修飾の情報のインスタンス
    formula : Formula
        オリゴマーの組成式のインスタンス
    deltams : float
        組成不明の質量差
    formula_str : str
        オリゴマーの組成式の文字列
    mol_w : mol_w
        オリゴの分子量
    ex_ms : float
        オリゴの精密質量
    maxnum_activeproton ：int
        活性プロトンの数。Na塩の分子量を求める時、m/zを計算する時に使用する.
    extinction : int
        モル吸光係数
    nmol_od : float
        nmol/OD, ODからnmolに変換する係数
    ug_od : float
        µg/OD, ODからµgに変換する係数.
    sodiumform : Oligomer
        FLP Na塩のインスタンス

    SeeAlso
    -------
    Compound

    """

    def __init__(
        self,
        constants: ConstantsForSequence,
        sequence: Sequence,
        modification_5: str or Modification = "H",
        modification_3: str or Modification = "H",
        name: str | None = None,
        note: str | None = None,
        deltams: float = 0.0,
        extinction: bool = False,
        sodiumform: bool = False,
        **kwargs,
    ):
        """初期化メソッド.

        Parameters
        ----------
        constants : ConstantsForSequence
            定数をまとめたインスタンス.
        sequence : Sequence
            配列情報のインスタンス
        modification_5 : Modification, optional The default is "H".
            5'末端の修飾情報のインスタンス
        modification_3 : Modification, optional The default is "H".
            5'末端の修飾情報のインスタンス
        name : str, optional The default is None.
            オリゴの名称.
        note : str, optional The default is None.
            配列の補足情報
        deltams : float, optional The default is 0.0.
            組成不明の質量差
        extinction : bool, optional The default is False.
             モル吸光係数関係の値を計算するか。
        sodiumform : bool, optional The default is False.
             Na塩のオリゴ情報を計算するか。
        **kwargs : TYPE
            表から読み取ったときに余計な列があるとエラーになるのでこれで受けておく.

        """
        self.constants = constants
        self.sequence = sequence
        if modification_5 == "H":
            self.modification_5 = constants.dict_modifications5["04-00-00"]
        elif isinstance(modification_5, str):
            self.modification_5 = constants.dict_modifications5[modification_5]
        else:
            self.modification_5 = modification_5
        if modification_3 == "H":
            self.modification_3 = constants.dict_modifications3["05-00-00"]
        elif isinstance(modification_3, str):
            self.modification_3 = constants.dict_modifications3[modification_3]
        else:
            self.modification_3 = modification_3
        self.formula = (
            self.sequence.get_formula() + self.modification_5.formula + self.modification_3.formula
        )
        super().__init__(
            name=name, note=note, deltams=deltams, formula=self.formula, constants=constants
        )
        self.maxnum_activeproton = self.count_activeproton()
        # モル吸光係数
        if extinction:
            self.set_extinction()
        else:
            self.extinction = 0
            self.nmol_od = 0
            self.ug_od = 0
        # Na塩の計算
        if sodiumform:
            self.sodiumform = Oligomer.get_sodiumform(self)
        else:
            self.sodiumform = None

    def count_activeproton(self) -> int:
        """活性プロトンの数をカウントする（m/zの価数maxやNa塩への変換に使用する値）."""
        activeproton = (
            self.sequence.length - 1 + self.modification_5.salt + self.modification_3.salt
        )
        return activeproton

    def set_extinction(self):
        """ectinction系のデータをまとめて更新."""
        self.extinction = self.sequence.get_extinction()
        self.extinction += self.modification_5.extinction
        self.extinction += self.modification_3.extinction
        self.nmol_od = 1 / self.extinction * 10**6
        self.ug_od = self.mol_w / self.extinction * 10**3

    def get_sodiumform(self):
        """FLPのNa塩のインスタンスを得る."""
        sodiumform = copy.deepcopy(self)
        sodiumform.formula.H -= self.maxnum_activeproton
        sodiumform.formula.Na += self.maxnum_activeproton
        sodiumform.mol_w = sodiumform.formula.get_ms(mol_w=True, constants=self.constants)
        sodiumform.ex_ms = sodiumform.formula.get_ms(mol_w=False, constants=self.constants)
        sodiumform.set_extinction()
        sodiumform.formula_str = str(sodiumform.formula)
        return sodiumform

    def expand_for_summary(self, mz_data=False, mol_w=True) -> dict:
        """FLP SummarySheet出力用に情報を抜きだす関数."""
        dict_data = {
            "Name": self.name,
            "Note": self.note,
            "Length": self.sequence.length,
            "Modification5": self.modification_5.name,
            "Sequence": self.sequence.sequence_ns,
            "Modification3": self.modification_3.name,
            "Formula": self.formula_str,
            "Ex.MS": self.ex_ms,
            "Mol.W": self.mol_w,
        }
        if self.extinction != 0:
            dict_data.update(
                {
                    "Extinction": self.extinction,
                    "nmol/OD": self.nmol_od,
                    "µmol/OD": self.ug_od,
                }
            )
        if isinstance(self.sodiumform, self.__class__):
            dict_data.update(
                {
                    "Formula(Na)": self.sodiumform.formula_str,
                    "Ex.MS(Na)": self.sodiumform.ex_ms,
                    "Mol.W(Na)": self.sodiumform.mol_w,
                }
            )
        if mz_data:
            dict_mz = self.make_mz_dict(num_z=self.maxnum_activeproton, mol_w=mol_w)
            dict_data.update(dict_mz)
        else:
            pass
        return dict_data


if __name__ == "__main__":
    # ForDEBUG
    const = ConstantsForSequence(
        filename_constant=FILENAME_CONSTANT,
        filename_sequence=FILENAME_SEQUENCE,
    )
    AA = Sequence.create_fm_sequence_ns("C*A*G*T*C*A*G*T", const)
    BB = Oligomer(constants=const, sequence=AA, extinction=True)
    CC = BB.expand_for_summary(mol_w=False, mz_data=True)
