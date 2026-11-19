# xpsloss-starter

公開理論に基づく、XPSの誘電損失計算の独立実装に向けた着手版です。
作成日: 2026-11-19。

**現時点で動くのは、次の4つです。**
- **振動子で表した損失関数**
- **その因果的な複素誘電応答**
- **均質・等方な無限媒質のバルク非弾性散乱**
- **法線出射・固定深さでの表面・正孔・真空中の項（Simonsen 1997 の公表式。ギャップのない材料のみ）**

**深さ平均、斜出射、ギャップのある材料での表面・正孔の項、スペクトルは未実装です。QUEELS-XPSの完全再現版ではありません。分散のある損失関数では、QUEELS-XPS v2.1 の出力とは一致しません（[理由](docs/QUEELS_DISPERSION_TEST.md)）。**

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
| 文献由来の固定入力（ELFデータ集のSi・Au・Al、manualのSi）と来歴の記録 | 同梱。値の物理的妥当性は未検証 |
| 法線出射・固定深さの表面・真空中の損失・静止正孔・干渉（ギャップなし、[Simonsen 1997 式(13)–(17)](docs/M1_EQUATION_MAP.md)） | 実装済み・テストあり・別担当がレビュー。原著の図3・5・8を再現 |
| 同上の QUEELS-XPS v2.1 との照合 | 分散0の入力では一致。分散ありでは不一致。参照ソフトの実部の処方が公表された式と違うことを、事前登録試験で確認（[QUEELS_DISPERSION_TEST.md](docs/QUEELS_DISPERSION_TEST.md)） |
| 斜出射、ギャップのある材料での表面・正孔の項 | 未実装 |
| 深さ平均と自己無撞着な平均自由行程 | 未実装 |
| 多重散乱、スペクトル畳み込み、装置関数 | 未実装 |
| Quanty接続、CoO物性値、フィッティング | 未実装 |
| 元のQUEELS-XPSとの数値一致 | 主張しない（上記） |

ELFは Im[-1/epsilon]、DIIMFPは単位経路長・単位損失エネルギー当たりの非弾性散乱確率密度です。ELF単体を、測定されるXPSの損失確率と同一視しないでください。

## Windows PowerShellでの開始

Python 3.10以上を想定します。展開した `xpsloss-starter` フォルダーを開きます。

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

## AIエージェントとの共同開発

最初に `AGENTS.md`、`docs/PROJECT_PLAN.md`、`docs/EQUATIONS.md` を渡してください。`CLAUDE.md`にも同じ資料への案内があります。

次に依頼する仕事は `docs/NEXT_TASK.md` に書いてあります。バルクの実装を足場として、まず1997年原著の数式と計算手順の対応を確定します。ソフトの機能が増えたことと、物理モデルの妥当性が検証されたことは、別々に記録します。

## 検証と資料

- `docs/VALIDATION.md`: この配布物で実行した検証と未検証事項。
- `docs/EQUATIONS.md`: 数式、単位、検証の根拠。
- `docs/M1_EQUATION_MAP.md`: 1997年原著の式(13)–(17)の転記、記号と単位、項ごとの実装予定、未確定の点。
- `docs/M2_FIRST_CASE.md`: M2の最初のケース（Al、1000 eV、法線出射）、参照ソフトとの比較の許容幅（比較の前に固定）。
- `docs/REVIEW_CHECKLISTS.md`: 実装前・実装後・参照比較のレビュー項目。
- `docs/QUEELS_DISPERSION_TEST.md`: QUEELS-XPS v2.1 の分散の扱いを判別した事前登録試験と、その結果。
- `docs/DIELECTRIC_DESIGN.md`: 複素誘電応答の方式比較、符号・幅の規約、ギャップ端と極の扱い。
- `docs/BENCHMARK_INPUTS.md`: 固定入力の出典、利用条件、光学和則の監査。
- `docs/PROJECT_PLAN.md`: 完全モデルとQuanty接続へ進むための設計。
- `docs/REFERENCES.md`: 一次資料と用途。
- `NOTICES.md`: 独立実装、同梱データの帰属表示、公開時の整理。

Windows 11（Python 3.12）とLinux（Python 3.12、M0のみ）で確認済みです。macOSは未確認です。詳細は `docs/VALIDATION.md` を参照してください。
