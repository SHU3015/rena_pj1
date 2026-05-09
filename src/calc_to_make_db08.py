import datetime
from functools import wraps
from pathlib import Path
import pandas as pd
from peptide import Peptide, ConstantsForPeptide


DIR_NAME = Path(__file__).resolve().parent
FILENAME = DIR_NAME.parent / "データベース/08_peptide.xlsx"
LIGAND_DATA_PATH = DIR_NAME.parent / "データベース/data/ligand.csv"


def stop_watch(func):
    """工程の経過時間出力用Wrapper."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        """時間出力部."""
        print(f"\n>> Start   : {func.__name__:^16} |")
        start = datetime.datetime.now()
        result = func(*args, **kwargs)
        deltatime = datetime.datetime.now() - start
        minitues = int(deltatime.total_seconds() / 60)
        sec = deltatime.seconds - 60 * int(deltatime.total_seconds() / 60)
        msec = f"{deltatime.microseconds:06d}"[:-4]
        print(
            f">> End     : {func.__name__:^16} | ",
            f"RAP   : {minitues:02d}′{sec:02d}″{msec:}",
        )
        return result

    return wrapper


class Sheet:
    """Peptideの計算用データをのせたクラス."""

    def __init__(self, filename, sheetname):
        self.sheet = pd.read_excel(filename, sheetname).fillna(
            {"deltams": 0, "cyclo": 0, "n_terminal": "H", "c_terminal": "OH"}
        )
        self.constants = ConstantsForPeptide(filename)

    @stop_watch
    def calc_ms(self):
        """Peptideを計算する."""
        self.sheet["obj"] = self.sheet.apply(
            lambda x: Peptide.create_fm_1letter(
                name=x["name"],
                note=x["note"],
                deltams=x["deltams"],
                n_terminal=x["n_terminal"],
                id_sequence=x["id"],
                c_terminal=x["c_terminal"],
                sequence=x["sequence"],
                cyclo=x["cyclo"],
                constants=self.constants,
            ),
            axis=1,
        )

    def output_data(self):
        """Peptideの計算結果を出力する."""
        self.result = self.sheet.apply(
            lambda x: x["obj"].expand_data_for_peptide(), axis=1, result_type="expand"
        )
        self.result.to_csv(
            f"{DIR_NAME.parent}/データベース/data/08_peptide.csv",
            encoding="utf-8-sig",
            index=False,
        )

    def output_liganddata(self):
        """Peptideの計算結果を出力する."""
        self.result = self.sheet.apply(
            lambda x: x["obj"].expand_data_for_ligand(), axis=1, result_type="expand"
        )
        self.result.to_csv(LIGAND_DATA_PATH, encoding="utf-8-sig")


if __name__ == "__main__":
    sheet = Sheet(FILENAME, "sequence")
    sheet.calc_ms()
    sheet.output_data()
    sheet = Sheet(FILENAME, "ligand")
    sheet.calc_ms()
    sheet.output_liganddata()
