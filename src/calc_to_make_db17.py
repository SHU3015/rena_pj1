"""06_rawdatabaseから07_baselinesequenceを作成するためのプログラム."""

from pathlib import Path
import pandas as pd
import os
from calc_to_make_db02 import Sequence, SequenceForRena, ConstantsForSequence


DIR_NAME = Path(__file__).resolve().parent

for_const = f"{DIR_NAME.parent}/データベース/02_ss_oligonucleotide.xlsx"
main_sheet = f"{DIR_NAME.parent}/データベース/17_sequence_baseline_forPB.xlsx"

const = ConstantsForSequence(filename_sequence=for_const)


def make_baseline(line):
    if line["type_modification"] == "SO":
        baseline = Sequence.create_fm_sequence_ns(
            line["baseline"], const).get_complementary_baseline(
                rna=False, reverse=False
                )
    else:
        baseline = line["baseline"]
    oligo = SequenceForRena.create_fm_baseline_as_gapmer(
        sequence_baseline=baseline,
        sugar_wing=line["sugar_wing"],
        sugar_gap=line["sugar_gap"],
        linker_wing=line["link_wing"],
        linker_gap=line["link_gap"],
        wing=(line["wing5"], line["wing3"]),
        constants=const
        )
    return oligo


df_sheet = pd.read_excel(main_sheet, sheet_name="converter").fillna("")
df_sheet["obj"] = df_sheet.apply(make_baseline, axis=1)
df_sheet["sequence_ns"] = df_sheet["obj"].apply(lambda x: x.sequence_ns)
df_sheet["sequence_gd"] = df_sheet["obj"].apply(lambda x: x.get_seq_gd())

df_for_output = df_sheet[[
    "ID", 'Name', 'note', 'id_sequence', 'baselinename', 'baseline',
    'len', 'type_modification', 'sugar_wing', 'sugar_gap', 'link_wing',
    'link_gap', 'wing5', 'wing3', 'sequence_ns', 'sequence_gd', '手動調整？'
    ]]

df_for_output.to_csv(
    f"{DIR_NAME.parent}/データベース/data/calcd_DB17_baselinesequence.csv",
    encoding="utf-8-sig", index=False)
os.startfile(rf"{DIR_NAME.parent}/データベース/17_sequence_baseline_forPB.xlsx")
