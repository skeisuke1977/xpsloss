# xpsloss

公開理論に基づくXPSの誘電損失計算を、Python APIとCLIで再現するためのパッケージです。v0.1.0 の対象は M0、M1、M2-1 です。[English README](README.md)

**現時点で動くのは、次の4つです。**
- **振動子で表した損失関数**
- **その因果的な複素誘電応答**
- **均質・等方な無限媒質のバルク非弾性散乱**
- **法線出射・固定深さでの表面・正孔・真空中の項（Simonsen 1997 の公表式。ギャップのない材料のみ）**

**深さ平均、斜出射、ギャップのある材料での表面・正孔の項、スペクトルは未実装です。QUEELS-XPSの完全再現版ではありません。試験した分散のある入力では、QUEELS-XPS v2.1 の出力と一致しませんでした（[追試の条件と結果](docs/QUEELS_DISPERSION_TEST.md)）。**

このパッケージはQUEELS/QUASESの公式製品ではなく、その実行ファイル・DLL・ソースコードを利用していません。公開されている数式から新規に実装しています。

## 現在の範囲

| 機能 | 状態 |
| --- | --- |
| JSON入力の検証 | 実装済み |
| Drude-Lindhard型振動子からELFを計算 | 実装済み |
| 運動量分散を含むバルクDIIMFP | 実装済み |
| CSV、JSONによる出力、計算条件の記録 | 実装済み |
| 解析解・和則・収束との照合 | テストあり |
| 振動子ELFからの因果的な複素誘電応答（ギャップ時はKKで実部を再構成） | 実装済み・テストあり・別担当がレビュー |
| 表形式ELFからのKK再構成（参照・診断用） | 実装済み・テストあり・別担当がレビュー |
| 光学和則・一次モーメントの診断（入力は再規格化しない） | 実装済み・テストあり |
| 文献由来の固定入力（ELFデータ集のSi・Au・Al）と来歴の記録 | 同梱。値の物理的妥当性は未検証 |
| 法線出射・固定深さの表面・真空中の損失・静止正孔・干渉（ギャップなし、[Simonsen 1997 式(13)–(17)](docs/M1_EQUATION_MAP.md)） | 実装済み・テストあり・別担当のエージェントがレビュー。原著の図3・5・8は、別担当のエージェントが原著の式(20)を使った独立なコードで、図の読取り誤差の範囲で再現 |
| 同上の QUEELS-XPS v2.1 との照合 | 試験した分散0の合成入力では、全成分が最大 `total` の 0.1–0.4 % 以内で一致。試験した分散ありの入力では不一致。試験した条件の参照ソフトの出力が、実部の分子に分散前のエネルギーを使う処方と一致することを、事前登録試験で確認（[QUEELS_DISPERSION_TEST.md](docs/QUEELS_DISPERSION_TEST.md)） |
| 斜出射、ギャップのある材料での表面・正孔の項 | 未実装 |
| 深さ平均と自己無撞着な平均自由行程 | 未実装 |
| 多重散乱、スペクトル畳み込み、装置関数 | 未実装 |
| Quanty接続、CoO物性値、フィッティング | 未実装 |
| 元のQUEELS-XPSとの数値一致 | 主張しない（上記） |

ELFは Im[-1/epsilon]、DIIMFPは単位経路長・単位損失エネルギー当たりの非弾性散乱確率密度です。ELF単体を、測定されるXPSの損失確率と同一視しないでください。

## Windows PowerShellでの開始

Python 3.10以上を使います。展開した `xpsloss` のルートフォルダーを開きます。

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m xpsloss capabilities --json
.\.venv\Scripts\python.exe -m xpsloss validate examples/synthetic_bulk.json
.\.venv\Scripts\python.exe -m xpsloss elf examples/synthetic_bulk.json --out runs/elf-demo
.\.venv\Scripts\python.exe -m xpsloss dielectric examples/synthetic_bulk.json --q-inv-A 0.5 --out runs/dielectric-demo
.\.venv\Scripts\python.exe -m xpsloss bulk examples/synthetic_bulk.json --out runs/bulk-demo
.\.venv\Scripts\python.exe -m xpsloss xps-fixed-depth examples/synthetic_bulk.json --depth-A 10 --epsilon-momentum k_total --out runs/xps-fixed-depth-demo
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

`xps-fixed-depth` は、表面法線からの出射角0°、生成深さ `--depth-A` で、正孔ありの有効断面積と、その「正孔なし」「正孔（単位電荷）」成分、項ごとの値を出力します。`--epsilon-momentum` は必須です。`k_total`（全運動量で ε を評価）を選ぶと、深い極限でバルクの式に一致します。`k_parallel` は原著の閉形式の導出どおりの評価です。

文献由来の固定入力と実行例は [benchmarks/README.md](benchmarks/README.md)、出典・帰属表示・和則の検算は [docs/BENCHMARK_INPUTS.md](docs/BENCHMARK_INPUTS.md)、複素誘電応答の設計と規約は [docs/DIELECTRIC_DESIGN.md](docs/DIELECTRIC_DESIGN.md) にあります。

仮想環境のactivateは不要です。通常のインストールには依存パッケージ取得のため通信が必要ですが、計算処理そのものは外部APIを使いません。依存関係を事前に導入すればオフラインで計算できます。

## macOS / Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m xpsloss capabilities --json
.venv/bin/python -m xpsloss bulk examples/synthetic_bulk.json --out runs/bulk-demo
.venv/bin/python -m unittest discover -s tests -v
```

**`examples/` の例題は合成材料です。Al、CoO、Siなどの実測・文献データとして利用しないでください。** 文献由来の入力は `benchmarks/inputs/` にだけ置き、出典ごとに分けています。既存の出力ディレクトリには上書きしません。再計算するときは新しい出力名を指定します。

計算結果の `run.json` には、コマンドラインで指定した表記のままの入力パス `input_file` と、入力ファイルのバイト列の SHA-256 `input_sha256` が入ります。

## ライセンスと引用

Keisuke Sato がコードを MIT、自作文書を CC BY 4.0 で提供します。文書中の第三者の引用は各権利者の条件に従います。第三者データ（`benchmarks/sources/**` と `benchmarks/inputs/*elfdb*`）は両方の適用範囲から除外され、N. Pauly、F. Yubero、S. Tougaard による *ELF dielectric functions for various materials*, v1.0（2022）の CC BY 4.0 が適用されます。[MITライセンス](LICENSE)・[帰属表示](NOTICES.md)・[ELFデータ集](https://doi.org/10.5281/zenodo.6024064)・[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

ソフトウェアの引用には Keisuke Sato、xpsloss、バージョン 0.1.0 を記載してください。ソフトウェア DOI は公開時に追加します。データ集を使用した場合はその出典も引用してください。[引用情報](CITATION.cff)

## 検証と資料

- `docs/VALIDATION.md`: この配布物で実行した検証と未検証事項。
- `docs/EQUATIONS.md`: 数式、単位、検証の根拠。
- `docs/M1_EQUATION_MAP.md`: 1997年原著の式(13)–(17)の転記、記号と単位、項ごとの実装との対応、未確定の点。
- `docs/M2_FIRST_CASE.md`: M2の最初のケース（Al、1000 eV、法線出射）、参照ソフトとの比較の許容幅（比較の前に固定）。
- `docs/REVIEW_CHECKLISTS.md`: 実装前・実装後・参照比較のレビュー項目。
- `docs/QUEELS_DISPERSION_TEST.md`: 分散のある入力で、QUEELS-XPS v2.1 の出力が公表式と候補式のどちらに一致するかを照合した事前登録試験。内部の実装は未確認。
- `docs/ROADMAP.md`: 公開までのロードマップと、初版についての決定。
- `docs/NEXT_TASK.md`: 公開候補の検証と、その後のM2-2の作業条件。
- `docs/DIELECTRIC_DESIGN.md`: 複素誘電応答の方式比較、符号・幅の規約、ギャップ端と極の扱い。
- `docs/BENCHMARK_INPUTS.md`: 固定入力の出典、利用条件、光学和則の監査。
- `docs/PROJECT_PLAN.md`: 完全モデルとQuanty接続へ進むための設計。
- `docs/REFERENCES.md`: 一次資料と用途。
- `NOTICES.md`: ライセンスの適用範囲、同梱データの帰属表示、参照ソフトとの関係。

Windows 11（Python 3.12）とLinux（Python 3.12、M0のみ）で確認済みです。macOS、Python 3.10 と 3.13 は未確認です。詳細は `docs/VALIDATION.md` を参照してください。
