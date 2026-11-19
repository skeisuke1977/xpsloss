# ベンチマークの固定入力

`inputs/` に出典ごとの入力記録4件と、データ集のSi・Au・Alを使ったバルク実行設定3件を収録しています。入力記録には原典の位置、利用条件、材料パラメータを保存しています。マニュアルのSiは相対強度で、絶対強度への換算手順が未確認のため実行設定はありません。[出典、帰属表示、和則の検算](../docs/BENCHMARK_INPUTS.md)

`sources/zenodo6024064_v1.0/` は Pauly、Yubero、Tougaard のデータ集 Version 1.0 の原本です。入力記録のSHA-256は改行を含む原本のバイト列に対応します。[データ集（DOI 10.5281/zenodo.6024064）](https://doi.org/10.5281/zenodo.6024064)、[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## CLIで確認する

リポジトリのルートで、NumPyとSciPyを導入したPythonを使います。作業ツリーのコードを確実に使うため、PowerShellでは `PYTHONPATH` を設定します。

```powershell
$env:PYTHONPATH = 'src'
python -m xpsloss validate benchmarks/inputs/si_elfdb_v1.0.bulk.json
python -m xpsloss validate benchmarks/inputs/au_elfdb_v1.0.bulk.json
python -m xpsloss validate benchmarks/inputs/al_elfdb_v1.0.bulk.json
python -m xpsloss bulk benchmarks/inputs/al_elfdb_v1.0.bulk.json --out runs/al-elfdb-bulk
```

`*.bulk.json` は、任意の動作確認条件として初期運動エネルギー1000 eV、損失格子1–100 eV（1 eV刻み）を指定します。バルクモデルは屈折率と内殻励起補正の `ac`、`bc` を使いません。QUEELS-XPSの数値出力との比較は別途行います。

実行結果の `run.json` の `config.provenance` には `source_id`、`reference_input`、`run_conditions_source`、`strength_basis` が記録されます。任意の `scope_note` と `raw_file_sha256` も保存されます。実行設定に来歴を付ける場合、強度区分は `absolute_eV2` が必須です。設定の解析では申告された値の区分、型と書式を検証します。参照先ファイルを開かないため、別ディレクトリに移した設定も読み込めます。同梱の入力記録と原本のハッシュの照合はテストで行います。

## 原本を読み取る

```python
from xpsloss.elfsource import read_elf

source = read_elf("benchmarks/sources/zenodo6024064_v1.0/ELF.txt", material="Al")
al, = source.materials
```

`read_elf` は数値を `Decimal` として保持し、実行時の浮動小数点精度で有限か確認します。ISO-8859-1とBOM付きUTF-8を読み取り、LF・CRLF・CRを行区切りとして扱います。文献一覧がある場合は、読み取る材料の文献番号が一覧にあるか検証します。材料名を指定すると、その区画を検証し、ファイル全体のSHA-256を返します。材料名を省略すると全区画を検証します。原本 `ELF.txt` のPPY区画で `@` が欠落しているため、全件読込では415行目のエラーが返ります。

## 参照ソフトとの比較（M2-1）

`compare_queels_fixed_depth.py` は、利用者が手元で計算した QUEELS-XPS の出力を読み、[M2_FIRST_CASE.md](../docs/M2_FIRST_CASE.md) §4.2 で事前に固定した基準で M2-1 の結果と比べます。出力するのは差の指標だけです。参照ソフトの出力は同梱せず、公開もしません（[NOTICES.md](../NOTICES.md)）。

```powershell
$env:PYTHONPATH = 'src'
python benchmarks/compare_queels_fixed_depth.py --reference-dir <QUEELS-XPSの出力フォルダー> --out work/comparisons/m2_1_al_1000_0.json
```

2026-11-20 の比較の結果は、分散のあるAlで不合格でした。原因の特定は [QUEELS_DISPERSION_TEST.md](../docs/QUEELS_DISPERSION_TEST.md) にあります。
