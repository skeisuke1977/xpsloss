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

実行結果の `run.json` の `config.provenance` に、`source_id`、`reference_input`、`run_conditions_source`、任意の `raw_file_sha256` が記録されます。`provenance` を指定する場合は最初の3項目が空でない文字列で、ハッシュを添える場合は小文字の16進数64桁です。原本との照合はテストで行います。

## 原本を読み取る

```python
from xpsloss.elfsource import read_elf

source = read_elf("benchmarks/sources/zenodo6024064_v1.0/ELF.txt", material="Al")
al, = source.materials
```

`read_elf` は数値を `Decimal` として保持します。材料名を指定すると、その区画を検証し、ファイル全体のSHA-256を返します。材料名を省略すると全区画を検証します。原本 `ELF.txt` のPPY区画で `@` が欠落しているため、全件読込では415行目のエラーが返ります。
